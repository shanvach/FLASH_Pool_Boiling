FROM akashdhruv/flash:4.0

MAINTAINER adhruv 

ARG flash_object
ARG flash_run_files
ARG flash_path

COPY $flash_path /home/object

COPY docker/Makefile.h /home/object/$flash_object

WORKDIR /home/object/$flash_object

RUN make && \
    mkdir /home/run/ && \ 
    cp flash.par /home/run/ && \
    cp flash4 /home/run/

WORKDIR /home/run

CMD ["sh", "-c", "./flash4"]
