# svn $Id: Module.mk 29 2007-04-23 19:23:26Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2013 BFM System Team (bfm_st@lists.cmcc.it)
# Copyright (c) 2002-2007 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := Build/BFM

local_lib  := libBFM.a
local_src  := $(wildcard $(local_sub)/*.?90)

$(eval $(call make-library,$(local_lib),$(local_src)))

$(eval $(compile-rules))
