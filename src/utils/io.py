import json
import os
import sys

class Data:
    @staticmethod
    def saveJSON(path: str, fileName: str, data: dict, readable=True) -> None:
        """Save a dictionary as a JSON file."""
        out_path = f"{path}/"
        try:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            with open(out_path + fileName, "w") as f:

                if readable:
                    json_str = json.dumps(data, indent=2)
                    import re

                    json_str = re.sub(r"\[\s+(\d+)\s+\]", r"[\1]", json_str)
                    # Keep simple two-element arrays on one line
                    json_str = re.sub(
                        r"\[\s+(\d+),\s+(\d+)\s+\]", r"[\1, \2]", json_str
                    )
                    f.write(json_str)
                # Format JSON more compactly while keeping it readable
                # json_str = json.dumps(data, indent=2)
                else:

                    json.dump(data, f, separators=(",", ":"))
                # Keep simple single-element arrays on one line
                # import re
                # json_str = re.sub(r'\[\s+(\d+)\s+\]', r'[\1]', json_str)
                # Keep simple two-element arrays on one line
                # json_str = re.sub(r'\[\s+(\d+),\s+(\d+)\s+\]', r'[\1, \2]', json_str)
                # f.write(json_str)
        except Exception as e:
            Logger.cprint(f"Failed to save JSON to {path}: {e}", Colors.FAIL)

    @staticmethod
    def loadJSON(path: str) -> dict:
        """Load a dictionary from a JSON file."""
        try:
            with open(path, "r") as f:
                return json.load(f)
        except Exception as e:
            Logger.cprint(f"Failed to load JSON from {path}: {e}", Colors.FAIL)
            return {}
        Logger.cprint("-" * 100, color)
