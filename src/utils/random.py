import uuid


def get_next_submission_id() -> str:
    return uuid.uuid4().hex


def get_random_uuid_hex():
    return uuid.uuid4().hex
