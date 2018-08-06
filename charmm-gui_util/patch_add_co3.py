#!/usr/bin/env python
"""
Patch CHARMM .inp files so we can add CO3 to it as anion
"""

__author__ = 'Michael King <michael.king@uni-konstanz.de'


import os

###############################################################################
# co3.crd
###############################################################################

co3_crd="""\
* CHARMM coordinates generated from VMD
         4  EXT
         1         1  CO3       C1              0.0000000000        0.0000000000        0.0000000000  HETB      11              0.0000000000
         2         1  CO3       O1             -1.0840001106        0.6259994507        0.3839988708  HETB      11              0.0000000000
         3         1  CO3       O2              0.8449997902       -0.4629993439        0.8799991608  HETB      11              0.0000000000
         4         1  CO3       O3              0.2920000553       -0.0340003967       -1.2470016479  HETB      11              0.0000000000
"""

with open("co3.crd", "w") as fp:
  fp.write(co3_crd)

###############################################################################
# step2.2_ions.inp
###############################################################################

INP = "step2.2_ions.inp"

print("Process: {}".format(INP))

inp_file = open(INP, 'r').read()
if not os.path.exists(INP+".bck"):
    with open(INP+".bck", "w") as fp:
        fp.write(inp_file)
else:
    print("Use: {}".format(INP+".bck"))
    inp_file = open(INP+".bck", 'r').read()

#=====================================================#
# define ncarbonate
find='set negval   = 1\n'
sub = find + "set carbonateval = 2\n"
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# calc ncarbonate
find = 'calc nneg = @nneg + int ( @conc * 6.021 * 0.0001 * @volumn ) * @posval\n'
sub = find + '''\
calc ncarbonate = int( @nneg / @carbonateval )
calc nneg = @nneg - @ncarbonate * @carbonateval
'''
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# calc diff
find = 'calc diff   = int ( ?cgtot + @npos*@posval - @nneg*@negval )'
sub  = 'calc diff   = int ( ?cgtot + @npos*@posval - @nneg*@negval - @ncarbonate * @carbonateval )'
assert len(inp_file.replace(find,sub)) != len(inp_file), find
inp_file=inp_file.replace(find,sub)

#=====================================================#
# count nneg
find = 'if diff .gt. 0 calc nneg = @nneg + 1\n'
sub = '''\
   if diff .gt. 0 then 
    calc mod = @diff - 2 * int( @diff / 2)
    if mod .gt. 0 then
      calc nneg = @nneg + 1
    else
      calc ncarbonate = @ncarbonate + 1
    endif
   endif
'''
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# Generate CO3
find='! Randomly place the ions\n!\n'
sub=find+'''
!Generate CO3
if @ncarbonate .gt. 0 then
  read sequence CO3 @ncarbonate
  generate CO3 first none last none setup warn
  open read card unit 10 name co3.crd
  read coor card unit 10 append
  
  set i 2
  label loop
  if i .le. @ncarbonate then
    coor dupl select segid CO3 .and. resid 1 end select segid CO3 .and. resid @i end
    increase i by 1
    goto loop
  endif
endif
'''
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# calc nion
find = 'calc nion = @npos + @nneg\n'
sub  = '''\
calc nion = @npos + @nneg + @ncarbonate
calc nwoneg = @npos + @ncarbonate
'''
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# calc xpos ypos zpos
find= '''\
calc xpos = @A / 2.0
calc ypos = @B / 2.0
calc zpos = @C / 2.0'''
sub = '''\
calc xpos = 0 ! @A / 2.0
calc ypos = 0 ! @B / 2.0
calc zpos = 0 ! @C / 2.0
'''
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)


#=====================================================#
# fix selection
find = 'cons fix sele .not. ( segid CAL .or. segid CLA ) end'
sub  = 'cons fix sele .not. ( segid CAL .or. segid CLA .or. segid CO3 ) end'
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)




#=====================================================#
# decide CO3
find = 'set ion = CLA\n       calc j  = @j - @npos\n'
sub  = '''\
if j .gt. @nwoneg then 
         set ion = CLA
         calc j  = @j - @nwoneg
       else
         set ion = CO3
         calc j  = @j - @npos
       endif
'''
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)
#=====================================================#
# coor set (TEST)

find = 'coor set xdir @xpos  ydir @ypos  zdir @zpos select target end\n'
sub  = '!' + find + '''\
    coor copy comp
    if ion .eq. CO3 then
      calc phi = 360 * ?random
      coor rota phi @phi xdir ?random ydir ?random zdir ?random xcen 0 ycen 0 zcen 0 select target end
    endif
    coor tran comp xdir -@xpos ydir -@ypos zdir -@zpos select target end
'''
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# dist calculation
find = 'coor dist cut'
sub  = 'coor dist comp cut'
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# coor set (SAVE)
find = "coor set xdir @xsave  ydir @ysave  zdir @zsave select target end\n       goto doinit\n"
sub = "!" + find + """\
    else
       coor tran xdir -@xpos ydir -@ypos zdir -@zpos  select target end 
       update
"""
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# !Image centering by residue
find = '.or. resname CLA end'
sub  = '.or. resname CLA .or. resname CO3 end'
assert len(inp_file.replace(find,sub)) != len(inp_file), find
inp_file=inp_file.replace(find,sub)

#=====================================================#
# set fix
find = 'cons fix sele .not. ( segid CAL .or. segi CLA ) end'
sub  = 'cons fix sele .not. ( segid CAL .or. segid CLA .or. segid CO3 ) end'
assert len(inp_file.replace(find,sub)) != len(inp_file), find
inp_file=inp_file.replace(find,sub)

#=====================================================#
# set fix
find = 'delete atom sele .not. ( segid CAL .or. segid CLA ) end'
sub  = 'delete atom sele .not. ( segid CAL .or. segid CLA .or. segid CO3 ) end'
assert len(inp_file.replace(find,sub)) != len(inp_file), find
inp_file=inp_file.replace(find,sub)


#=====================================================#
# set OUTPUT
# number 
find = '* set nneg = @nneg'
sub  = '* set ncarbonate = @ncarbonate ! Number of carbonate ions\n' + find
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

# id set
find = '* set negid = CLA'
sub  = '* set carbid = CO3\n' + find
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

# sequence
find = '* read sequence CLA'
sub  = '* read sequence CO3 @ncarbonate !Generate CO3\n* generate CO3 first none last none setup\n' + find
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

with open(INP, 'w') as fp:
    fp.write(inp_file)
    
###############################################################################
# step2_solvator.inp
###############################################################################


INP = "step2_solvator.inp"
print("Process: {}".format(INP))

inp_file = open(INP, 'r').read()
if not os.path.exists(INP+".bck"):
    with open(INP+".bck", "w") as fp:
        fp.write(inp_file)
else:
    print("Use: {}".format(INP+".bck"))
    inp_file = open(INP+".bck", 'r').read()

#=====================================================#
# define ncarbonate
find='stream step2.2_ions.prm\n'
sub = find + """
if ncarbonate .gt. 0 then
   read sequence @carbid @ncarbonate
   generate CO3 first none last none setup warn
endif
"""
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# calc nion
find = 'calc nion = @npos + @nneg\n'
sub  = 'calc nion = @npos + @nneg + @ncarbonate\n'
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

#=====================================================#
# set OUTPUT
# number 
find = '* set nneg = @nneg'
sub  = '* set ncarbonate = @ncarbonate ! Number of carbonate ions\n' + find
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

# id set
find = '* set negid = @negid'
sub  = '* set carbid = @carbid\n' + find
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

# define
find = '* if npos .gt. 0 then\n'
sub = """\
* if ncarbonate .gt. 0 then
*  read sequence @carbid @ncarbonate
*  generate CO3 first none last none setup warn
* endif
""" + find
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)

with open(INP, 'w') as fp:
    fp.write(inp_file)


###############################################################################
# step3_pbcsetup.inp
###############################################################################

INP = "step3_pbcsetup.inp"
print("Process: {}".format(INP))

inp_file = open(INP, 'r').read()
if not os.path.exists(INP+".bck"):
    with open(INP+".bck", "w") as fp:
        fp.write(inp_file)
else:
    print("Use: {}".format(INP+".bck"))
    inp_file = open(INP+".bck", 'r').read()


#=====================================================#
# segid carbid
# number 
find = 'segid @posid .or. segid @negid'
sub  = 'segid @posid .or. segid @negid .or. segid @carbid'
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub)

#=====================================================#
# set OUTPUT
# id set
find = '* set negid = @negid'
sub  = '* set carbid = @carbid\n' + find
assert len(inp_file.replace(find,sub,1)) != len(inp_file), find
inp_file=inp_file.replace(find,sub,1)


with open(INP, 'w') as fp:
    fp.write(inp_file)