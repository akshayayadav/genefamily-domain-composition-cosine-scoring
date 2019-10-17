FROM ubuntu:16.04
RUN apt-get update && \
    apt-get install -y wget build-essential make perl hmmer && \
    cd /root/ && \
    wget "http://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/OldPfamScan/PfamScan1.5/PfamScan.tar.gz" -O- \
    | tar --no-same-owner -xz && \
    cpan -fi Moose  && \
    cpan -fi IPC::Run && \
    cpan -fi Bio::SeqFeature::Generic && \
    cd  PfamScan/ && \
    mv Bio/ /usr/share/perl/5.22/ && \
    mv pfam_scan.pl /usr/local/bin/ && \
    cd /root && \
    rm -rf PfamScan* && \
    apt-get purge -y build-essential make && \
    apt-get autoremove -y --purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && pip3 install --upgrade pip \
  && apt-get autoremove \
  && apt-get clean

RUN pip3 install 'snakemake==5.5.4'
RUN pip3 install pandas
RUN pip3 install matplotlib
RUN pip3 install sklearn

ADD scripts /usr/local/bin
RUN chmod +x /usr/local/bin/*
