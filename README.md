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

## Deploy into production

1. Add Gitlab CI Variable `AWS_ECR_REPOSITORY_NAME` for AWS ECR initialization on CI/CD.
2. Add `helm/values.yaml` to set up service configs for Kubernetes deployment (e.g [values.yaml](https://git.ndfs.kz/ndfs/infra/tools/deployers/kubernetes-deployer/-/blob/main/deployer/helm/values/example.yaml)). 
3. Add `.gitlab-ci.yml` to set up CI/CD. For more information about CI/CD jobs, visit [gitlab-ci templates](https://git.ndfs.kz/ndfs/templates/gitlab-ci).
```yaml
include:
  - project: 'ndfs/templates/gitlab-ci'
    ref: main
    file: 'service-fastapi.gitlab-ci.yml'

deployer-kubernetes:
  variables:
    DEPLOYER_HELM_NAME: '<unqiue project name>'
    DEPLOYER_HELM_VALUES: '<path to helm values, e.g helm/values.yaml>'
```
