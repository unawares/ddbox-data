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
2. Add `helm/values.yaml` to set up service configs for Kubernetes deployment.
```yaml
global:
  namespace: <set-me>
  defaultImage: <set-me>
  secrets: {
    # <secret-name>:
    #   SECRET_KEY: <secret-key>
  }
  environments: {
    # ENV_NAME: ENV_VALUE
  }

app:
  name: <set-me>
  # !!! Uncomment: put configs here if you want to use !!!
  # configs:
  #   base.json: |
  #     {
  #       "env": "dev"
  #     }

  # You can remove ingress if there is no need of outside access
  ingress:
    hosts:
    - paths:
      - prefixPath: /
        servicePort: <set-me>  # set the service port
    annotations: {}

  service:
    # !!! Uncomment configs if you want to set it !!!
    # configs:
    # - mountPath: /configs/base.json
    #   subPath: base.json
    replicas: 1
    ports:
    - portName: main
      servicePort: <set-me>  # set the service port (used by ingress)
      containerPort: <set-me>  # set the port mapping
    # !!! Uncomment: change CMD command off a docker image if it is necessary
    # command: <set-me>
    resources:
      requests:
        memory: 64Mi
        cpu: 0.04
      limits:
        memory: 256Mi
        cpu: 0.16
```
3. Add `.gitlab-ci.yml` to set up CI/CD. For more information about CI/CD jobs, visit [gitlab-ci templates](https://git.ndfs.kz/ndfs/templates/gitlab-ci).
```yaml
stages:
  - init
  - test
  - scan
  - build
  - deploy

include:
  - project: 'ndfs/templates/gitlab-ci'
    ref: main
    file: 'aws-ecr-init.gitlab-ci.yml'
  - project: 'ndfs/templates/gitlab-ci'
    ref: main
    file: 'tester-pytest.gitlab-ci.yml'
  - project: 'ndfs/templates/gitlab-ci'
    ref: main
    file: 'scanner-semgrep.gitlab-ci.yml'
  - project: 'ndfs/templates/gitlab-ci'
    ref: main
    file: 'scanner-safety.gitlab-ci.yml'
  - project: 'ndfs/templates/gitlab-ci'
    ref: main
    file: 'builder-fastapi.gitlab-ci.yml'
  - project: 'ndfs/templates/gitlab-ci'
    ref: main
    file: 'deployer-kubernetes.gitlab-ci.yml'
  - project: 'ndfs/templates/gitlab-ci'
    ref: main
    file: 'optimizations.gitlab-ci.yml'

.base: &base
  only:
    refs:
      - main
  tags:
    - ndfs

aws-ecr-init:
  <<: *base
  stage: init

tester-pytest:
  <<: *base
  stage: test
  only:
    refs:
      - main
      - merge_requests

scanner-semgrep:
  <<: *base
  stage: scan
  only:
    refs:
      - main
      - merge_requests

scanner-safety:
  <<: *base
  stage: scan
  only:
    refs:
      - main
      - merge_requests

builder-fastapi:
  <<: *base
  stage: build

deployer-kubernetes:
  <<: *base
  stage: deploy
  variables:
    DEPLOYER_HELM_NAME: '<unqiue project name>'
    DEPLOYER_HELM_VALUES: '<path to helm values, e.g helm/values.yaml>'
```
