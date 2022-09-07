help: ## display this help message
	@echo "Please use \`make <target>' where <target> is one of"
	@grep '^[a-zA-Z]' ${MAKEFILE_LIST} | sort | awk -F ':.*?## ' 'NF==2 {printf "\033[36m  %-25s\033[0m %s\n", $$1, $$2}'

install: ## install
	poetry install

dev.start: install ## start
	poetry run uvicorn src.app.main:app --host 0.0.0.0 --port 8000 --reload

requirements.txt: install ## generate requirements.txt
	poetry export -f requirements.txt --output requirements.txt

# Flags

.PHONY: *
