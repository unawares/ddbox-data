from typing import List, Dict
import requests


def submit_moses(smiles_list: List[str]):
    api_url = 'http://localhost:8000/submission/moses'
    response = requests.post(api_url, json={
        'smiles_list': smiles_list
    })

    if response.status_code != 200:
        raise Exception(response.text)
    
    return response.json()['submission_id']


def submit_docking(smiles_list: List[str], receptors: dict):
    api_url = 'http://localhost:8000/submission/docking'

    response = requests.post(api_url, json={
        'smiles_list': smiles_list,
        'receptor_ids': [receptor['id'] for receptor in receptors],
        'centers': [receptor['center'] for receptor in receptors],
        'sizes': [receptor['size'] for receptor in receptors],
    })

    if response.status_code != 200:
        raise Exception(response.text)
    
    return response.json()['submission_id']


def submission_moses_status(submission_id: str):
    api_url = 'http://localhost:8000/submission/result/moses'
    response = requests.get(api_url, params={
        'submission_id': submission_id,
    })

    if response.status_code != 200:
        raise Exception(response.text)

    return response.json()


def submission_docking_status(submission_id: str):
    api_url = 'http://localhost:8000/submission/result/docking'
    response = requests.get(api_url, params={
        'submission_id': submission_id,
    })

    if response.status_code != 200:
        raise Exception(response.text)

    return response.json()
