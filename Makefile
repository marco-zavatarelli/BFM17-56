##-----------------------------------------------------------------------
## Makefile for making a new distribution release
##-----------------------------------------------------------------------

BFM_RELEASE="4.0"
SYSTEM_RELEASE=""

all: distclean 

distclean:
	make -C src/ $@
	@$(RM)  ./bin/* 

current_release: distclean
	@scripts/release.sh $(BFM_RELEASE) $(SYSTEM_RELEASE) current

release: distclean
	@scripts/release.sh $(BFM_RELEASE) $(SYSTEM_RELEASE) 

##-----------------------------------------------------------------------
##
##-----------------------------------------------------------------------
