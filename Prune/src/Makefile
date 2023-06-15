CC=gcc
CXX=g++
PREFIX=../

all: htslib
	$(CXX) -O3 --std=c++11 -o ALLHiC_prune -I. -Ihtslib-1.17 ALLHiC_prune.cpp Prune.cpp htslib-1.17/libhts.a -lz -llzma -lbz2 -lpthread -lcurl

htslib: untar
	$(MAKE) -C htslib-1.17

untar:
	@if [ ! -d "htslib-1.17" ]; then echo "Untar"; tar jxvf htslib-1.17.tar.bz2; else echo "Skip"; fi

clean:
	rm -rf ALLHiC_prune htslib-1.17

install: all
	@if [ ! -d $(PREFIX) ]; then mkdir -p $(PREFIX); fi
	cp ALLHiC_prune $(PREFIX)
