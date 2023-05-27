import json
from typing import *


def json_to_pretty_str(json_data: Mapping[str, Any]):
    json_formatted_str = json.dumps(json_data, indent=2)
    return json_formatted_str
