From python:latest

LABEL maintainer="James W Gillespie"
LABEL maintainer-email="gillejw@auburn.edu"

RUN mkdir /usr/src/pdbt
WORKDIR /usr/src/pdbt

ENV PYTHONUNBUFFERED 1
