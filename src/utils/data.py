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

    def formatted(self) -> str:
        return format(self.value, self.fmt) if self.fmt else str(self.value)

    def __str__(self) -> str:
        return f"{self.label}: {self.formatted()}"

    def __repr__(self) -> str:
        return str(self)
