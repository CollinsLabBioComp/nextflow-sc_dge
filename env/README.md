
# Setting up environment

We use [conda](https://docs.conda.io/en/latest/) and [renv](https://rstudio.github.io/renv/index.html) for managing necessary packages. **Note**: we recommend using [mamba](https://mamba.readthedocs.io/en/latest/) instead of base conda to speed up download times. If mamba is not available, conda should still work.

## 1. Setting up conda environment

Install the required packages via conda/mamba:
```bash
# The repo directory.
REPO_MODULE="${HOME}/repo/path/to/this/pipeline"

# Install environment using Mamba. Replace 'mamba' with 'conda' if mamba not available.
mamba env create --name sc_diff_expr --file ${REPO_MODULE}/env/environment.yml

# Activate the new Conda environment.
source activate sc_diff_expr

# To update environment file:
#conda env export --no-builds | grep -v prefix | grep -v name > environment.yml
```

## 2. Setting up R environment

We install base R through the conda environment, but for R package management we use `renv`. `renv` is installed through the conda file, so we just need to load the rest of the R packages.

```bash
cd ${REPO_MODULE}/env
R
```

Now, we load the necessary packages using `renv`
```R
# renv::restore() will load packages detailed in renv.lock file
renv::restore(lockfile = 'renv.lock')

# It will ask whether or not to activate a profile. If you choose 'No', a profile will be activated in home directory. If you press 'Yes', the profile will activate in current directory.
# You can specify location of profile (location of cache to download packages) like:
renv::restore(lockfile = 'renv.lock', project='/path/to/cache/dir')
```

As an FYI, to build a `renv` profile:
```R
# initialize renv
renv::init()

# update renv after installing any packages
renv::settings$snapshot.type("simple")
renv::snapshot()

# add renv.lock to git. Note: Do not add .Rprofile and renv/activate.R to git. This will force R to use the renv directory and not conda meaning the packages will not work in the pipeline.
```

If `renv` does not work, the packages are listed in [rpackage_details.txt](rpackage_details.txt) and can be manually installed with:
```r
install.packages("<package>")
```

## 3. [OPTIONAL] Docker

Alternatively, we use the Docker image of the environment [https://hub.docker.com/repository/docker/letaylor/sc_nf_diffexpression](https://hub.docker.com/repository/docker/letaylor/sc_nf_diffexpression).
