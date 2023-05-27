from typing import *


class SuccessResponce:

    def __init__(self, response) -> None:
        self.response = response

    def to_response(self):
        return self.response
