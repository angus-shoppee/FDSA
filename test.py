
from splice import SamFlag


flag = SamFlag(2211)

print("FLAG:", flag)
print("Not primary alignment:", SamFlag.NOT_PRIMARY_ALIGNMENT in flag)
