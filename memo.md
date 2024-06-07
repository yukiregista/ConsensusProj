Instructions for auto-generating the document

1. At the project root, run `sphinx-apidoc -f -M -o docs/source src/Consensus`. 
This will create rst files at docs/source from Consensus module.

2. Go to `docs` folder, then run `make html`.
  - You might need to `make clean` to remove build cache to apply changes (eps. if they are made to conf.py)

3. For publishing with github pages, the html files needs to be under `docs/`. For that, we can use the following command:
    `sphinx-build -b html ./docs/source ./docs `