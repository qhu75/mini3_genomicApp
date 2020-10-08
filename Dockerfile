FROM virtualstaticvoid/heroku-docker-r:shiny
ENV PORT=8080
RUN apt-get install -y libssl-dev && apt-get clean
CMD ["/usr/bin/R", "--no-save", "--gui-none", "-f", "/app/run.R"]
