from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Callable
from .data import ResultData

class Colors:
    # Bright colors
    HEADER = "\033[95m"
    BRIGHT_MAGENTA = "\033[95m"
    BRIGHT_CYAN = "\033[96m"
    BRIGHT_BLUE = "\033[94m"
    BRIGHT_GREEN = "\033[92m"
    BRIGHT_YELLOW = "\033[93m"
    BRIGHT_RED = "\033[91m"
    BRIGHT_WHITE = "\033[97m"

    # Standard colors
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    TURQUOISE = "\033[36m"
    WHITE = "\033[37m"

    # Named aliases
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    SUCCESS = "\033[92m"
    INFO = "\033[94m"

    # Styles
    BOLD = "\033[1m"
    DIM = "\033[2m"
    ITALIC = "\033[3m"
    UNDERLINE = "\033[4m"

    # Reset
    ENDC = "\033[0m"


class Logger:

    HEADLINE = Colors.BRIGHT_CYAN + Colors.BOLD
    INFO = Colors.BRIGHT_BLUE
    SUCCESS = Colors.BRIGHT_GREEN
    WARNING = Colors.BRIGHT_YELLOW
    ERROR = Colors.BRIGHT_RED
    NOTICE = Colors.BRIGHT_MAGENTA
    NORMAL = Colors.ENDC

    @staticmethod
    def cprint(message: str, color: str = "", bold: bool = False) -> None:
        style = color
        if bold:
            style += Colors.BOLD
        print(f"{style}{message}{Colors.ENDC}")

    @staticmethod
    def header(message: str, color: str = "") -> None:
        Logger.cprint("-" * 100, color)
        Logger.cprint(message, color, True)
        Logger.cprint("-" * 100, color)

    @staticmethod
    def print_results(
        rows: list[list[ResultData]],
        color_fn: Callable[[list[ResultData]], str] | None = None,
    ) -> None:
        """Print a list of result rows as a right-justified table.

        rows     : each inner list is one row; labels from the first row
                   become the column headers
        color_fn : optional callable(row) -> Colors.XYZ string applied per row;
                   defaults to Colors.BOLD when omitted
        """
        if not rows:
            return

        visible_indices = [i for i, rd in enumerate(rows[0]) if rd.show]
        if not visible_indices:
            return

        headers = ["|  " +  rows[0][i].label + " " for i in visible_indices]
        widths = [len(h) for h in headers]

        # Pre-format every cell and compute column widths
        all_formatted: list[list[str]] = []
        for row in rows:
            fmted = []
            for col, i in enumerate(visible_indices):
                rd = row[i]
                s = rd.formatted()
                fmted.append(s)
                display_len = rd.width if rd.width > 0 else len(s)
                widths[col] = max(widths[col], display_len)
            all_formatted.append(fmted)

        header_str = " ".join(h.rjust(w) for h, w in zip(headers, widths))
        print(header_str)
        print("-" * len(header_str))

        for row, fmted in zip(rows, all_formatted):
            clr = Colors.FAIL if any(rd.is_flagged() for rd in row) else (color_fn(row) if color_fn else Colors.BOLD)
            parts = []
            for col, (i, s) in enumerate(zip(visible_indices, fmted)):
                rd = row[i]
                if rd.width > 0:
                    # Value already has its own padding/color — just add left-pad to align
                    parts.append(" " * max(0, widths[col] - rd.width) + s)
                else:
                    parts.append(s.rjust(widths[col]))
            print(f"{clr}{' '.join(parts)}{Colors.ENDC}")
