# Running GillespieSSA2 in a Docker container

[rcannood/gillespiessa2](https://hub.docker.com/r/rcannood/gillespiessa2) contains all necessary packages to run GillespieSSA2 from start to finish.

## Running the container
To run the container, you can use the following command.

```sh
docker run --rm -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true -v `pwd`:/home/rstudio/workdir rcannood/gillespiessa2
```

<!-- fedora users currently need to run:
docker run --rm --ulimit="nofile=4096" -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true -v `pwd`:/home/rstudio/workdir rcannood/gillespiessa2
-->

Keep this window open, and open up a browser and go to [127.0.0.1:8787](127.0.0.1:8787). Open up the file `vignettes/an_introduction.Rmd` or any of the other vignettes for a small example on how to run SSA simulations with GillespieSSA2.

The command can be dissected as follows.

```sh
docker run \

  # remove container after use
  --rm \
  
  # specify which port rstudio server uses
  -p 127.0.0.1:8787:8787 \
  
  # disable authentication because I'm lazy
  -e DISABLE_AUTH=true \
  
  # mount the current working directory to the rstudio home folder
  # so you can save results in the 'workdir' folder for later use
  -v `pwd`:/home/rstudio/workdir/ \
  
  # specify which container to run
  rcannood/gillespiessa2
```

## Update the container

If a newer version of the container has been released, you can update it by running the following command.
```sh
docker pull rcannood/gillespiessa2
```


## Building the container

To rebuild this docker container from scratch, run the following command with GillespieSSA2 (from the GillespieSSA2 folder, not the docker folder).

```sh
docker build -t rcannood/gillespiessa2 -f docker/Dockerfile --build-arg GITHUB_PAT=$GITHUB_PAT .
```

GITHUB_PAT should be an environment variable corresponding to the Personal Access Token created by following [this tutorial](https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token).


