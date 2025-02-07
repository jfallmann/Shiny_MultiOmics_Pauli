FROM rocker/shiny-verse:4.4.1
ARG DEBIAN_FRONTEND=noninteractive
ENV SHINY_LOG_STDERR=1
ENV TZ=UTC

RUN apt-get update \     
    && apt-get clean \
    && apt-get autoremove -y \
    && apt-get autoclean -y 
    # && apt-get install -y tzdata wget libxml2-dev \
    # && libfontconfig1-dev&& ln -s /usr/share/zoneinfo/Europe/Vienna /etc/localtime 

# Application specific packages
RUN mkdir -p /packages/R \
    && mkdir -p /srv/data \
    && mkdir -p /srv/docs
COPY R /packages/R
COPY app/data /srv/data
COPY app/docs /srv/docs

ENV R_REMOTES_UPGRADE="never"
RUN	R -e "install.packages('renv', repos = 'https://cloud.r-project.org/')" \
    && R -e "renv::restore('/packages/R')"

RUN chown shiny:shiny /var/lib/shiny-server \ 
    && rm -rf /srv/shiny-server/* \
    && chown shiny:shiny /srv/shiny-server

ADD --chown=shiny:shiny ./app/mitozebra.R /srv/shiny-server/mitozebra.R
ADD --chown=shiny:shiny ./app/server.R /srv/shiny-server/server.R
ADD --chown=shiny:shiny ./app/ui.R /srv/shiny-server/ui.R

EXPOSE 3838
USER shiny

CMD R -e 'shiny::runApp("/srv/shiny-server/mitozebra.R", port = 3838, host = "0.0.0.0")'
