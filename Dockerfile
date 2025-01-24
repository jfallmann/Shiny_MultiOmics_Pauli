FROM registry.hub.docker.com/library/ubuntu:24.04
ARG DEBIAN_FRONTEND=noninteractive
ENV SHINY_LOG_STDERR=1

RUN apt-get update && \
    ln -s /usr/share/zoneinfo/Europe/Vienna /etc/localtime && \
    apt-get install -y tzdata r-base gdebi-core wget libxml2-dev libfontconfig1-dev && \
    R -e "install.packages('shiny', repos='https://cran.rstudio.com/')" && \
    wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.22.1017-amd64.deb && \
    gdebi -n shiny-server-1.5.22.1017-amd64.deb && \
    rm shiny-server-1.5.22.1017-amd64.deb

# Application specific packages
RUN mkdir -p /packages/R
COPY R /packages/R

ENV R_REMOTES_UPGRADE="never"
RUN	R -e "source('/packages/R/setup_R_env.R')"

RUN chown shiny:shiny /var/lib/shiny-server && \
    rm -rf /srv/shiny-server/* && \
    chown shiny:shiny /srv/shiny-server

ADD --chown=shiny:shiny ./entrypoint.sh /entrypoint.sh
ADD --chown=shiny:shiny ./app /srv/shiny-server
RUN chmod +x /entrypoint.sh

EXPOSE 3838
USER shiny

# CMD R -e 'shiny::runApp("/srv/shiny-server/server.R", port = 3838, host = "0.0.0.0")'
ENTRYPOINT ["/entrypoint.sh"]