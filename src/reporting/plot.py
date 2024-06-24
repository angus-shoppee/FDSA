# Feature-Directed Splice Analysis
# Copyright (C) 2024 Angus Shoppee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import base64
from typing import Dict, List, Tuple, Union
from io import BytesIO
from statistics import mean
from pandas import DataFrame
from matplotlib import pyplot as plt
from matplotlib import patches
from matplotlib.path import Path

from reporting.process_results import QuantifiedSpliceJunctionLocus, FdsaResult


def find_location_within_exon(
    nucleotide: int,
    exon_positions: List[Tuple[int, int]]
) -> Union[Tuple[None, None], Tuple[int, int]]:
    """Takes a nucleotide position and returns its exon number and distance from exon start"""

    for i, position in enumerate(exon_positions):
        if position[0] <= nucleotide <= position[1]:
            return i + 1, nucleotide - position[0]
    return None, None  # nucleotide is not located within an exon


def find_flanking_exon_numbers(
    nucleotide: int,
    exon_positions: List[Tuple[int, int]],
    forward_stranded: bool
) -> Union[Tuple[None, None], Tuple[int, int]]:
    n_exons = len(exon_positions)
    for i, position in enumerate(exon_positions):
        if i + 1 == n_exons:
            continue
        next_position = exon_positions[i + 1]
        if forward_stranded:
            if position[1] < nucleotide < next_position[0]:
                return i + 1, i + 2
        else:
            if position[0] > nucleotide > next_position[1]:
                return i + 1, i + 2
    return None, None  # nucleotide is not located within an intron


def convert_global_exon_positions_to_local(
    exon_positions: List[Tuple[int, int]],
    forward_stranded: bool
) -> List[Tuple[int, int]]:
    n_exons = len(exon_positions)

    offset = exon_positions[0][0] if forward_stranded else exon_positions[-1][0]

    intron_distance: Dict[int: int] = {}

    if forward_stranded:

        for i, position in enumerate(exon_positions):

            # Distance is zero for first exon
            if i == 0:
                intron_distance[i + 1] = 0
                continue

            last_position = exon_positions[i - 1]

            intron_distance[i + 1] = position[0] - last_position[1] - 1
    else:

        for i, position in enumerate(exon_positions):

            # Distance is zero for final exon
            if i + 1 == n_exons:
                intron_distance[i + 1] = 0
                continue

            last_position = exon_positions[i + 1]

            intron_distance[i + 1] = position[0] - last_position[1] - 1

    scaled_exon_positions: List[Tuple[int, int]] = []
    _cumulative_intron_distance = 0
    # _exon_positions_range is coerced back to list to prevent consumption after iteration
    _exon_positions_range = range(n_exons) if forward_stranded else list(reversed(list(range(n_exons))))

    for i in _exon_positions_range:
        #     for i in range(n_exons)

        position = exon_positions[i]

        _cumulative_intron_distance += intron_distance[i + 1]

        scaled_position = (
            position[0] - offset - _cumulative_intron_distance,
            position[1] - offset - _cumulative_intron_distance
        )

        if forward_stranded:
            scaled_exon_positions.append(scaled_position)
        else:
            scaled_exon_positions.insert(0, scaled_position)

    return scaled_exon_positions


class OutOfBoundsError(ValueError):
    pass


class PlotJunction:
    data: QuantifiedSpliceJunctionLocus
    start_exon_nearest: Union[None, int]
    end_exon_nearest: Union[None, int]
    plot_start: int
    plot_end: int
    length: int
    exon_delta: int
    is_regular_splicing: bool
    highlight: bool

    def __init__(
        self,
        quantified_splice_junction_locus: QuantifiedSpliceJunctionLocus,
        global_exon_positions: List[Tuple[int, int]],
        plotting_exon_positions: List[Tuple[int, int]],
        intron_display_length: int,
        forward_stranded: bool,
        highlight: bool = False
    ) -> None:

        # TODO: Visualise junctions that originate outside of transcript bounds? (currently skipped)

        self.data = quantified_splice_junction_locus
        self.highlight = highlight

        start_exon, start_steps_from_exon_start = find_location_within_exon(
            quantified_splice_junction_locus.start - 1,  # Subtract 1 to reach prev. exon end
            global_exon_positions
        )

        if start_exon is None:
            _start_intron = find_flanking_exon_numbers(
                quantified_splice_junction_locus.start,
                global_exon_positions,
                forward_stranded
            )
            if forward_stranded:
                start_intron_left, start_intron_right = _start_intron
            else:
                start_intron_right, start_intron_left = _start_intron
            if start_intron_left is None or start_intron_right is None:
                raise OutOfBoundsError(
                    f"Junction location {quantified_splice_junction_locus.start} is out of bounds.\n" +
                    f"Exons: {global_exon_positions}"
                )
            _start_fraction_along_intron = (
                                               quantified_splice_junction_locus.start -
                                               global_exon_positions[start_intron_left - 1][1]
                                           ) / (global_exon_positions[start_intron_right - 1][0] -
                                                global_exon_positions[start_intron_left - 1][1])
            self.plot_start = \
                plotting_exon_positions[start_intron_left - 1][1] + \
                int(intron_display_length * _start_fraction_along_intron)
            # OLD / BUGGED (calculated position from start, instead of end, of exon to the left of intronic locus):
            # plotting_exon_positions[end_intron_left - 1][0] + \
        else:
            self.plot_start = \
                plotting_exon_positions[start_exon - 1][0] + start_steps_from_exon_start
            start_intron_left = None  # Keep the linter happy (start_intron_left is always defined)

        self.start_exon_nearest = start_exon if start_exon is not None else start_intron_left

        end_exon, end_steps_from_exon_start = find_location_within_exon(
            quantified_splice_junction_locus.end + 1,  # Add 1 to reach next exon start
            global_exon_positions
        )

        if end_exon is None:
            _end_intron = find_flanking_exon_numbers(
                quantified_splice_junction_locus.end,
                global_exon_positions,
                forward_stranded
            )
            if forward_stranded:
                end_intron_left, end_intron_right = _end_intron
            else:
                end_intron_right, end_intron_left = _end_intron
            if end_intron_left is None or end_intron_right is None:
                raise OutOfBoundsError(
                    f"Junction location {quantified_splice_junction_locus.end} is out of bounds.\n" +
                    f"Exons: {global_exon_positions}"
                )
            _end_fraction_along_intron = (
                                             quantified_splice_junction_locus.end -
                                             global_exon_positions[end_intron_left - 1][1]
                                         ) / (global_exon_positions[end_intron_right - 1][0] -
                                              global_exon_positions[end_intron_left - 1][1])
            self.plot_end = \
                plotting_exon_positions[end_intron_left - 1][1] + \
                int(intron_display_length * _end_fraction_along_intron)
            # OLD / BUGGED (calculated position from start, instead of end, of exon to the left of intronic locus):
            # plotting_exon_positions[end_intron_left - 1][0] + \
        else:
            self.plot_end = \
                plotting_exon_positions[end_exon - 1][0] + end_steps_from_exon_start
            end_intron_right = None  # Keep the linter happy (end_intron_right is always defined)

        self.end_exon_nearest = end_exon if end_exon is not None else end_intron_right

        # TODO: Raise a proper error here (though should never occur)
        assert all([self.start_exon_nearest, self.end_exon_nearest])  # Should never be None

        self.length = abs(self.plot_end - self.plot_start)

        _start_exon = start_exon if start_exon is not None else start_intron_left
        _end_exon = end_exon if end_exon is not None else end_intron_right
        self.exon_delta = max(1, abs(_end_exon - _start_exon))  # Give delta a floor value of 1

        self.is_regular_splicing = False
        if all([
            start_exon is not None,
            end_exon is not None,
            self.exon_delta == 1
        ]):
            if all([
                self.plot_start == plotting_exon_positions[start_exon - 1][1],
                self.plot_end == plotting_exon_positions[end_exon - 1][0]
            ]):
                self.is_regular_splicing = True


def plot_transcript(
    fdsa_result: FdsaResult,
    draw_junctions_with_min_n_occurrences: int = 1,
    show_main_title: bool = True,
    limit_to_samples: Union[None, List[str]] = None
) -> str:

    samples_to_plot = limit_to_samples if limit_to_samples is not None else fdsa_result.all.keys()

    n_samples = len(samples_to_plot)

    n_exons = len(fdsa_result.exon_positions)

    # TODO: Though an extreme edge case in the context of FDSA's intended purpose, this logic does not work for
    #       transcripts consisting of a single exon - should require strandedness as an explicit argument here
    #       (splice events should not be relevant for these, unless there is intra-exonic splicing)
    forward_stranded = True
    if len(fdsa_result.exon_positions) > 1:
        if fdsa_result.exon_positions[0][0] > fdsa_result.exon_positions[1][0]:
            forward_stranded = False

    avg_exon_length = mean([position[1] - position[0] for position in fdsa_result.exon_positions])

    scaled_exon_positions = convert_global_exon_positions_to_local(
        fdsa_result.exon_positions,
        forward_stranded
    )

    plotting_exon_positions = []
    plotting_spacer_scale = 0.67
    draw_intron_size = int(avg_exon_length * plotting_spacer_scale)
    _cumulative_plotting_spacer_between_exons = 0
    _exon_positions_range = range(n_exons) if forward_stranded else reversed(list(range(n_exons)))

    for i in _exon_positions_range:

        position_adjusted_for_plotting = (
            scaled_exon_positions[i][0] + _cumulative_plotting_spacer_between_exons,
            scaled_exon_positions[i][1] + _cumulative_plotting_spacer_between_exons
            # scaled_exon_positions[i][0] + int(_cumulative_plotting_spacer_between_exons * plotting_spacer_scale),
            # scaled_exon_positions[i][1] + int(_cumulative_plotting_spacer_between_exons * plotting_spacer_scale)
        )
        if forward_stranded:
            plotting_exon_positions.append(position_adjusted_for_plotting)
        else:
            plotting_exon_positions.insert(0, position_adjusted_for_plotting)

        _cumulative_plotting_spacer_between_exons += draw_intron_size
        # _cumulative_plotting_spacer_between_exons += avg_exon_length

    # ----------------------
    # DEFINE SHARED ELEMENTS

    rects = [
        [[(position[0], 1), position[1] - position[0], 1], "black", "none", 1]  # geometry, edgecolor, facecolor, alpha
        for position in plotting_exon_positions
    ]

    plot_feature_start_exon, plot_feature_start_steps = find_location_within_exon(
        fdsa_result.feature_region[0],
        fdsa_result.exon_positions
    )
    plot_feature_end_exon, plot_feature_end_steps = find_location_within_exon(
        fdsa_result.feature_region[1],
        fdsa_result.exon_positions
    )

    for exon_number in range(plot_feature_start_exon, plot_feature_end_exon + 1):
        # draw_rect_start = plotting_exon_positions[exon_number - 1][0] + plot_feature_start_steps \
        #     if exon_number == plot_feature_start_exon else plotting_exon_positions[exon_number - 1][0]
        # draw_rect_end = plotting_exon_positions[exon_number - 1][0] + plot_feature_end_steps \
        #     if exon_number == plot_feature_end_exon else plotting_exon_positions[exon_number - 1][1]
        draw_rect_start = plotting_exon_positions[exon_number - 1][0] + plot_feature_start_steps \
            if exon_number == plot_feature_start_exon else \
            (
                plotting_exon_positions[exon_number - 1][0]
                if forward_stranded
                else plotting_exon_positions[exon_number - 1][1]
            )
        draw_rect_end = plotting_exon_positions[exon_number - 1][0] + plot_feature_end_steps \
            if exon_number == plot_feature_end_exon else \
            (
                plotting_exon_positions[exon_number - 1][1]
                if forward_stranded
                else plotting_exon_positions[exon_number - 1][0]
            )

        rects.append([
            [(min(draw_rect_start, draw_rect_end), 1), abs(draw_rect_end - draw_rect_start), 1],  # geometry
            "none",  # edgecolor
            "red",  # facecolor
            0.5  # alpha
        ])

    xlim = [
        0 - avg_exon_length,
        (plotting_exon_positions[-1][1] if forward_stranded else plotting_exon_positions[0][1]) + avg_exon_length
    ]

    # -----------------------
    # DEFINE SUBPLOT ELEMENTS

    subplot_elements = {i: dict() for i in range(n_samples)}

    for plot_index, sample_name in enumerate(samples_to_plot):

        connections, curves, connection_annotations, curve_annotations = [], [], [], []

        connection_annotation_y_offset = 0.11
        curve_annotation_y_offset = 0.1
        curve_base_height = 0.5
        curve_height_multiplier = 0.75

        plot_junctions = []
        for junction in fdsa_result.all[sample_name]:
            _highlight = (junction in fdsa_result.overlapping[sample_name])
            # TODO: Only if not red
            if (
                junction.n >= draw_junctions_with_min_n_occurrences
            ) or (
                junction.n >= 1 and _highlight
            ):
                try:
                    plot_junctions.append(
                        PlotJunction(
                            junction,
                            fdsa_result.exon_positions,
                            plotting_exon_positions,
                            draw_intron_size,
                            forward_stranded,
                            highlight=_highlight
                        )
                    )
                except OutOfBoundsError:
                    # logger.warning("INFO: Skipped junction with out-of-bounds location")
                    continue

        # Order plot_junctions from the lowest length to highest
        plot_junctions.sort(key=lambda x: x.length)

        # Map each unique exon_delta value to a height value for plotting
        _all_exon_deltas_sorted = sorted(list({plot_junction.exon_delta for plot_junction in plot_junctions}))
        exon_delta_map_to_height_level = {delta: i + 1 for (i, delta) in enumerate(_all_exon_deltas_sorted)}

        # Keep track of similar junctions for application of a special height bonus to reduce overlapping
        encountered_any_spanning = {
            (plot_junction.start_exon_nearest, plot_junction.end_exon_nearest): 0
            for plot_junction in plot_junctions
        }

        # Keep track of highest curve to calculate ylim
        tallest_curve_reached = 0

        for plot_junction in plot_junctions:

            if plot_junction.is_regular_splicing:
                connections.append([
                    (plot_junction.plot_start, 1.5),
                    (plot_junction.plot_end, 1.5)
                ])
                connection_annotations.append([[
                    plot_junction.data.n,
                    (
                        mean((plot_junction.plot_start, plot_junction.plot_end)),
                        1.5 + connection_annotation_y_offset
                    )
                ], "black"])
            else:
                height_level = exon_delta_map_to_height_level[plot_junction.exon_delta] + \
                               encountered_any_spanning[
                                   (plot_junction.start_exon_nearest, plot_junction.end_exon_nearest)]
                curve_height = height_level * curve_base_height
                true_height = curve_height * curve_height_multiplier
                tallest_curve_reached = true_height if true_height > tallest_curve_reached else tallest_curve_reached
                curves.append([
                    [
                        (plot_junction.plot_start, 2),
                        (plot_junction.plot_start, 2 + curve_height),
                        (plot_junction.plot_end, 2 + curve_height),
                        (plot_junction.plot_end, 2)
                    ],  # geometry
                    "red" if plot_junction.highlight else "black",  # edgecolor
                    1.0 if plot_junction.highlight else 0.7  # alpha
                ])
                curve_annotations.append([[
                    plot_junction.data.n,
                    (
                        mean((plot_junction.plot_start, plot_junction.plot_end)),
                        2 + true_height + curve_annotation_y_offset
                    )
                ], "red" if plot_junction.highlight else "black"])

            if plot_junction.start_exon_nearest is not None and plot_junction.end_exon_nearest is not None:
                encountered_any_spanning[(plot_junction.start_exon_nearest, plot_junction.end_exon_nearest)] += 1

        exon_annotations = [
            [[
                f"E{i + 1}",
                (
                    mean((plotting_exon_positions[i][0], plotting_exon_positions[i][1])),
                    0.8 - (0.05 * tallest_curve_reached)
                )
            ], "black"]
            for i in range(n_exons)
        ]

        ylim = [0.5 - (0.03 * tallest_curve_reached), 2.5 + (1.05 * tallest_curve_reached)]

        subplot_elements[plot_index]["sample_name"] = sample_name
        subplot_elements[plot_index]["exon_annotations"] = exon_annotations
        subplot_elements[plot_index]["curve_annotations"] = curve_annotations
        subplot_elements[plot_index]["connection_annotations"] = connection_annotations
        subplot_elements[plot_index]["curves"] = curves
        subplot_elements[plot_index]["connections"] = connections
        subplot_elements[plot_index]["exon_annotations"] = exon_annotations
        subplot_elements[plot_index]["tallest_curve_reached"] = tallest_curve_reached
        subplot_elements[plot_index]["ylim"] = ylim

    # -------------
    # DRAW SUBPLOTS

    figure_width = 16

    fig, axs = plt.subplots(
        n_samples,
        figsize=(figure_width, sum([1 + subplot_elements[plot_index]["ylim"][1] for plot_index in range(n_samples)])),
        gridspec_kw={"height_ratios": [1 + subplot_elements[plot_index]["ylim"][1] for plot_index in range(n_samples)]}
    )

    for plot_index in range(n_samples):

        sample_name = subplot_elements[plot_index]["sample_name"]
        # exon_annotations = subplot_elements[plot_index]["exon_annotations"]
        curve_annotations = subplot_elements[plot_index]["curve_annotations"]
        connection_annotations = subplot_elements[plot_index]["connection_annotations"]
        curves = subplot_elements[plot_index]["curves"]
        connections = subplot_elements[plot_index]["connections"]
        exon_annotations = subplot_elements[plot_index]["exon_annotations"]
        tallest_curve_reached = subplot_elements[plot_index]["tallest_curve_reached"]
        ylim = subplot_elements[plot_index]["ylim"]

        ax = axs[plot_index]

        for rect in rects:
            ax.add_patch(
                patches.Rectangle(
                    *rect[0],
                    fill=True,
                    edgecolor=rect[1],
                    facecolor=rect[2],
                    alpha=rect[3],
                    lw=2
                )
            )
        for connection in connections:
            ax.add_patch(
                patches.PathPatch(
                    Path(connection, [Path.MOVETO, Path.LINETO]),
                    lw=2
                )
            )
        for curve in curves:
            ax.add_patch(
                patches.PathPatch(
                    Path(curve[0], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                    facecolor="none",
                    edgecolor=curve[1],
                    alpha=curve[2],
                    lw=2
                )
            )
        for annotation in exon_annotations + curve_annotations:
            ax.annotate(
                *annotation[0],
                ha="center",
                color=annotation[1],
                size=16
            )
        for annotation in connection_annotations:
            ax.text(
                annotation[0][1][0],
                annotation[0][1][1] + 0.015 * tallest_curve_reached,
                annotation[0][0],
                color=annotation[1],
                size=16,
                ha="center",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.5)
            )

        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        ax.axis("off")

        ax.set_title(sample_name, size=24)

    main_title = f"{fdsa_result.gene_name} - Feature {fdsa_result.feature_number} " + \
                 f"of {fdsa_result.total_features_in_transcript}"

    _first_ax_fig = axs[0].get_figure()

    if show_main_title:
        _first_ax_fig.suptitle(main_title, size=42)

    _first_ax_fig.tight_layout()
    # TODO: Fix main title height scaling across different plot heights
    _first_ax_fig.subplots_adjust(top=0.958)  # Works for 18 samples
    #     _first_ax_fig.subplots_adjust(top=0.962)  # Works for 14 samples

    # Write plot to bytes buffer
    svg_buffer = BytesIO()
    plt.savefig(
        svg_buffer,
        format="svg",
        bbox_inches="tight"
    )
    plt.close()

    # Return as base64 string
    return base64.b64encode(svg_buffer.getvalue()).decode("utf-8")


def plot_splice_rate(
    fdsa_result: FdsaResult,
    norm_gene_counts: DataFrame,
    sample_groups: Dict[str, List[str]],
    group_name_by_sample: Dict[str, str],
    color_by_group_name: Dict[str, str],
    shape_by_group_name: Dict[str, str],
    show_main_title: bool = True
) -> str:

    x, y, color, shape, labels = [], [], [], [], []
    for sample_name, frequency in fdsa_result.frequencies.items():
        y.append(frequency)
        x.append(norm_gene_counts.loc[fdsa_result.gene_id][sample_name])
        color.append(color_by_group_name[group_name_by_sample[sample_name]])
        shape.append(shape_by_group_name[group_name_by_sample[sample_name]])
        labels.append(group_name_by_sample[sample_name])

    fig, ax = plt.subplots(figsize=(8, 6))

    for i in range(len(x)):
        ax.scatter(
            x[i],
            y[i],
            c=color[i],
            marker=shape[i],
            label=labels[i],
            s=100
        )

    plt.xlabel("Gene expression (TMM norm. CPM)", fontsize=22)
    plt.ylabel("Splice event detection (%)", fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=18)

    handles, labels = ax.get_legend_handles_labels()
    handles_for_unique_labels, unique_labels = [], []
    # Legend entries will be ordered as defined in the config file in python versions 3.7 and above,
    # where dictionaries preserve order
    for handle, label in sorted(
        list(zip(handles, labels)),
        key=lambda v: list(sample_groups.keys()).index(v[1])
    ):
        if label not in unique_labels:
            unique_labels.append(label)
            handles_for_unique_labels.append(handle)

    plt.legend(
        handles_for_unique_labels,
        unique_labels,
        fontsize=18,
        frameon=False,
        # loc="right",
        # bbox_to_anchor=(1, 0.5)
        loc="center left",
        bbox_to_anchor=(1.02, 0.5)
    )

    legend = plt.gca().get_legend()

    if show_main_title:
        ax.set_title(
            f"{fdsa_result.gene_name} - Feature {fdsa_result.feature_number}" +
            f" of {fdsa_result.total_features_in_transcript}",
            size=28
        )

    # Write plot to bytes buffer
    svg_buffer = BytesIO()
    plt.savefig(
        svg_buffer,
        format="svg",
        bbox_extra_artists=(legend,),
        bbox_inches="tight"
    )
    plt.close()

    # Return as base64 string
    return base64.b64encode(svg_buffer.getvalue()).decode("utf-8")
