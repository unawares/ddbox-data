import uuid


def get_next_submission_id() -> str:
    return uuid.uuid4().hex
