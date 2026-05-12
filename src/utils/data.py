from dataclasses import dataclass
from typing import Any


@dataclass
class Data:
    label: str
    value: Any


@dataclass
class ResultData(Data):
    show: bool = True
    fmt: str = ""
    width: int = 0  # 0 = auto; set explicitly when value contains ANSI codes
    flag: str = ""  # "non-zero" | "zero" | "negative" — triggers error color on the row

    def formatted(self) -> str:
        if self.fmt == "factors":
            from nt.common import fmt_factored
            return fmt_factored(int(self.value))
        return format(self.value, self.fmt) if self.fmt else str(self.value)

    def is_flagged(self) -> bool:
        if self.flag == "non-zero":
            return self.value != 0
        if self.flag == "zero":
            return self.value == 0
        if self.flag == "negative":
            return self.value < 0
        return False

    def __str__(self) -> str:
        return f"{self.label}: {self.formatted()}"

    def __repr__(self) -> str:
        return str(self)
