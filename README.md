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
