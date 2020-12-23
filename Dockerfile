FROM debian:10-slim AS builder

RUN apt update && apt install -y --no-install-recommends \
  autoconf \
  automake \
  g++ \
  libtool \
  make \
  yaggo \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src

COPY . .

RUN autoreconf -fi \
  && ./configure \
  && make \
  && make install

# optionally run tests
#RUN make check || cat ./test-suite.log && exit 1

FROM debian:10-slim

RUN apt update && apt install -y --no-install-recommends \
  gnuplot \
  libgomp1 \
  perl \
  && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local /usr/local
RUN ldconfig
