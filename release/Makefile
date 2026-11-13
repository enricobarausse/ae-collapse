# Build the Einstein-aether collapse code.
#   make            -> ./collapse   (g++ -O2 -std=c++11, the flags used for all paper runs)
#   make warn       -> compile with -Wall -Wextra (no output expected)
#   make smoke      -> build, then a short run (N = 1601, A0 = 0.35, to t = 8; a few
#                      seconds) in the directory smoke_run/, as a sanity check
#   make clean
CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++11

collapse: collapse.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

warn: collapse.cpp
	$(CXX) $(CXXFLAGS) -Wall -Wextra -fsyntax-only $<

smoke: collapse
	mkdir -p smoke_run
	cd smoke_run && ../collapse --N=1601 --A0=0.35 --t_final=8 --times=2,4,8 > log.txt 2> err.txt
	@echo "last line of smoke_run/phi.dat (step t dt phi):"; tail -n 1 smoke_run/phi.dat; ls smoke_run

clean:
	rm -f collapse
	rm -rf smoke_run

.PHONY: warn smoke clean
