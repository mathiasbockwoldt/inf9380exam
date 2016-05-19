# Galaxy Haystack

FROM bgruening/galaxy-stable
MAINTAINER Mathias Bockwoldt mathias.bockwoldt@uit.no

WORKDIR /galaxy-central

RUN add-tool-shed --url 'http://testtoolshed.g2.bx.psu.edu/' --name 'Test Tool Shed'

# Install Visualisation
# RUN install-biojs msa

# Adding the tool definitions to the container
ADD my_tool_list.yml $GALAXY_ROOT/my_tool_list.yml

# Install Haystack
RUN install-tools $GALAXY_ROOT/my_tool_list.yml

# Mark folders as imported from the host.
VOLUME ["/export/", "/data/", "/var/lib/docker"]

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800

# Autostart script that is invoked during container start
CMD ["/usr/bin/startup"]
