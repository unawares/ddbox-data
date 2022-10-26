# FastAPI

## Initialization from [FastAPI template](https://git.ndfs.kz/ndfs/templates/stacks/fastapi)

Use the following commands to copy the template and set to a new repository:

```bash
git clone git@git.ndfs.kz:ndfs/templates/stacks/fastapi.git
cd fastapi
git remote rename origin template-origin
git remote add origin <repository-url>
git push --set-upstream origin main
```

Install [poetry](https://python-poetry.org/docs/):

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

Initialize the project:

```bash
make init
```

## Getting started

Use the help command to see the available commands in Makefile:
```bash
make help
```

Run the development server (the server is accessible at [http://localhost:8000](http://localhost:8000)):
```bash
make dev.start
```

Run the development server using docker compose (the server is accessible at [http://localhost:8000](http://localhost:8000)):
```bash
make docker.compose.up
```

## Deployment

1. Add Gitlab CI Variable `AWS_ECR_REPOSITORY_NAME` in `dev`/`prod` environment for AWS ECR initialization on CI/CD.
2. Add `helm/values.dev.yaml`/`helm/values.yaml` to set up service configs for Kubernetes deployment (e.g. [example.yaml](https://git.ndfs.kz/ndfs/infra/tools/deployers/kubernetes-deployer/-/blob/main/deployer/helm/values/example.yaml)). 
3. Add `.gitlab-ci.yml` with the following content to set up CI/CD jobs (visit [gitlab-ci template](https://git.ndfs.kz/ndfs/templates/gitlab-ci/-/blob/main/service-fastapi.gitlab-ci.yml) for more information about CI/CD jobs).
4. Generate `requirements.txt` using the command `make requirements.txt` and commit the file.
```yaml
include:
  - project: 'ndfs/templates/gitlab-ci'
    ref: main
    file: 'service-fastapi.gitlab-ci.yml'

dev-deployer-kubernetes:
  variables:
    DEPLOYER_HELM_NAME: '<unqiue project name>'
    DEPLOYER_HELM_VALUES: '<path to helm values, e.g. helm/values.dev.yaml>'

prod-deployer-kubernetes:
  variables:
    DEPLOYER_HELM_NAME: '<unqiue project name>'
    DEPLOYER_HELM_VALUES: '<path to helm values, e.g. helm/values.yaml>'
```

## Project structure

- `.flake8` - flake8 is a python tool that glues together pycodestyle, pyflakes, mccabe, and third-party plugins to check the style and quality of some python code.
- `.mypy.ini` - mypy is an optional static type checker for python that aims to combine the benefits of dynamic (or "duck") typing and static typing.
- `.gitignore` - a gitignore file specifies intentionally untracked files that git should ignore.
- `.isort.cfg` - a config file for isort. Isort is a Python utility / library to sort imports alphabetically, and automatically separated into sections and by type.
- `.pre-commit-config.yaml` - a config file for pre-commit which is installed after `make init`.
- `Dockerfile` - a dockerfile only for the local deployment. For the production deployment you should specify the `BUILDER_DOCKERFILE` environment for the CI/CD builder job.
- `docker-compose.yml` - a docker compose file only for the local deployment. For the production deployment it is not used.
- `Makefile` - a collection of useful scripts for development.
- `README.md` - a documentation about the project.
- `pyproject.toml` - a meta file for poetry.
- `poetry.lock` - a lock file for poetry.
- `requirements.txt` - a pip requirements file generated by the command `make requirements.txt`.
- `src/*` - source code.
- `devops/*` - infratructure staffs.
- `devops/commands/entrypoint.sh` - entrypoint script run on deployment.
- `devops/commands/start.sh` - starting script run on deployment.
