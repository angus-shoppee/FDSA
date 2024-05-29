
from typing import List, Union, Any
from functools import reduce


def remove_quotes(value: str) -> str:
    return value.replace("\"", "").replace("\'", "")


def flatten_nested_lists(nested_lists: List[List[Any]]) -> List[Any]:
    return reduce(lambda a, b: a + b, nested_lists)


def divide_or_default_zero(
    numerator: Union[int, float],
    denominator: Union[int, float]
) -> float:

    return 0.0 if not all((numerator, denominator)) else numerator / denominator
