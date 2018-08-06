#!/usr/bin/env python
"""
Patch CHARMM .inp files so we can add CO3 to it as anion
"""

__author__ = 'Michael King <michael.king@uni-konstanz.de'

import os
import re

# CO3 file
_STR_co3_crd = """\
* CHARMM coordinates generated from VMD
         4  EXT
         1         1  CO3       C1              0.0000000000        0.0000000000        0.0000000000  HETB      11              0.0000000000
         2         1  CO3       O1             -1.0840001106        0.6259994507        0.3839988708  HETB      11              0.0000000000
         3         1  CO3       O2              0.8449997902       -0.4629993439        0.8799991608  HETB      11              0.0000000000
         4         1  CO3       O3              0.2920000553       -0.0340003967       -1.2470016479  HETB      11              0.0000000000
"""


def write_co3(fname="co3.crd"):
    # type: (str) -> None
    """
    Function to write a crd file containing CO3

    Parameters
    ----------
    fname : str, optional
        filename
    """
    with open(fname, "w") as fp:
        fp.write(_STR_co3_crd)


def _get_subtype(astring):
    """
    Helper function to interpret the substype

    * `"a"` for append,
    * `"p"` for prepend
    * `"s"` for substitute (append = prepend = False)
    * `int` for the number of replacements

    Parameters
    ----------
    astring : str
        string to interpret.
        May contain `"a"` for append,
        `"p"` for prepend,
        `"s"` for substitute (append = prepend = False),
        `int` for the number of replacements
    Returns
    -------
    prepend : bool
        if it should be prepended
    append : bool
        if it should be appended
    count : int or None
        Counts how often it should be replaced
    """
    # check if we should append and/or prepend
    append = True if astring.find("a") != -1 else False
    prepend = True if astring.find("p") != -1 else False

    if astring.find("s") != -1:
        append = prepend = False

    # get the integer in the string
    str_match = re.search(r'\d+', astring)
    count = None if str_match is None else int(str_match.group())

    return prepend, append, count


def patch_file(fname, list_replacements, output=None, backup=True, verbose=True):
    # type: (str, list[list[str]], Union[str,None], bool, bool) -> None
    """
    Function to patch a file with a given list of replacements

    Parameters
    ----------
    fname : str
        filename of the patchfile
    list_replacements : List[str]
        List of replacement strings in the form of:

        `list(sublist_1, sublist_2, ..., sublist_n)`

        with `sublist_i` = `[sub_type, find_str, sub_str]`

        `sub_type` :
            * `"a"` for append,
            * `"p"` for prepend
            * `"s"` for substitute (append = prepend = False)
            * `int` for the number of replacements
    output : Union[str,None]
        output filename
    backup : bool, optional
        if a backup should be written or not `fname` -> `fname.bck`
    verbose : bool, optional
        Turn verbose on / off
    """

    if verbose:
        print("Process: {}".format(fname))

    # read the file
    inp_file = open(fname, 'r').read()

    # create a backup if required
    if backup:
        if not os.path.exists(fname + ".bck"):
            with open(fname + ".bck", "w") as fp:
                fp.write(inp_file)
        else:
            if verbose:
                print("Use: {}".format(fname + ".bck"))
            inp_file = open(fname + ".bck", 'r').read()

    # set output
    output = fname if output is None else output

    # iterate over the lis of replacements
    for line in list_replacements:
        prepend, append, count = _get_subtype(line[0])
        find = line[1]
        sub = line[2]

        if prepend:
            sub = sub + find
        if append:
            if prepend:
                sub = sub + line[3]
            else:
                sub = find + sub

        if count is None:
            replaced_file = inp_file.replace(find, sub)
        else:
            replaced_file = inp_file.replace(find, sub, count)
        assert len(replaced_file) != len(inp_file), "Could not find:\n{}".format(find)
        inp_file = replaced_file

    # write output
    with open(output, 'w') as fp:
        fp.write(inp_file)


if __name__ == '__main__':
    # write CO3 file
    write_co3(fname="co3.crd")

    # =====================================================#
    # step2.2_ions.inp

    INP = "step2.2_ions.inp"

    list_replacements = [
        [  # define ncarbonate
            'a1',
            'set negval   = 1\n',
            'set carbonateval = 2\n'
        ],
        [  # calc ncarbonate
            'a1',
            'calc nneg = @nneg + int ( @conc * 6.021 * 0.0001 * @volumn ) * @posval\n',
            '''\
calc ncarbonate = int( @nneg / @carbonateval )
calc nneg = @nneg - @ncarbonate * @carbonateval
'''
        ],
        [  # calc diff
            's',
            'calc diff   = int ( ?cgtot + @npos*@posval - @nneg*@negval )',
            'calc diff   = int ( ?cgtot + @npos*@posval - @nneg*@negval - @ncarbonate * @carbonateval )'
        ],
        [  # count nneg
            's1',
            'if diff .gt. 0 calc nneg = @nneg + 1\n',
            '''\
   if diff .gt. 0 then 
    calc mod = @diff - 2 * int( @diff / 2)
    if mod .gt. 0 then
      calc nneg = @nneg + 1
    else
      calc ncarbonate = @ncarbonate + 1
    endif
   endif
'''
        ],
        [  # Generate CO3
            'a1',
            '! Randomly place the ions\n!\n',
            '''
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
        ],
        [  # calc nion
            's1',
            'calc nion = @npos + @nneg\n',
            '''\
calc nion = @npos + @nneg + @ncarbonate
calc nwoneg = @npos + @ncarbonate
'''
        ],
        [  # calc xpos ypos zpos
            's1',
            '''\
calc xpos = @A / 2.0
calc ypos = @B / 2.0
calc zpos = @C / 2.0''',
            '''\
calc xpos = 0 ! @A / 2.0
calc ypos = 0 ! @B / 2.0
calc zpos = 0 ! @C / 2.0
'''
        ],
        [  # fix selection
            's1',
            'cons fix sele .not. ( segid CAL .or. segid CLA ) end',
            'cons fix sele .not. ( segid CAL .or. segid CLA .or. segid CO3 ) end'
        ],
        [  # decide CO3
            's1',
            'set ion = CLA\n       calc j  = @j - @npos\n',
            '''\
if j .gt. @nwoneg then 
         set ion = CLA
         calc j  = @j - @nwoneg
       else
         set ion = CO3
         calc j  = @j - @npos
       endif
'''
        ],
        [  # coor set (TEST)
            'ap1',
            'coor set xdir @xpos  ydir @ypos  zdir @zpos select target end\n',
            '!',
            '''\
    coor copy comp
    if ion .eq. CO3 then
      calc phi = 360 * ?random
      coor rota phi @phi xdir ?random ydir ?random zdir ?random xcen 0 ycen 0 zcen 0 select target end
    endif
    coor tran comp xdir -@xpos ydir -@ypos zdir -@zpos select target end
'''
        ],
        [  # dist calculation
            's1',
            'coor dist cut',
            'coor dist comp cut'
        ],
        [  # coor set (SAVE)
            'ap1',
            "coor set xdir @xsave  ydir @ysave  zdir @zsave select target end\n       goto doinit\n",
            '!',
            """\
    else
       coor tran xdir -@xpos ydir -@ypos zdir -@zpos  select target end 
       update
"""
        ],
        [  # !Image centering by residue
            's',
            '.or. resname CLA end',
            '.or. resname CLA .or. resname CO3 end'
        ],
        [  # set fix
            's',
            'cons fix sele .not. ( segid CAL .or. segi CLA ) end',
            'cons fix sele .not. ( segid CAL .or. segid CLA .or. segid CO3 ) end'
        ],
        [  # delete atom
            's',
            'delete atom sele .not. ( segid CAL .or. segid CLA ) end',
            'delete atom sele .not. ( segid CAL .or. segid CLA .or. segid CO3 ) end'
        ],
        [  # set OUTPUT (number)
            'p1',
            '* set nneg = @nneg',
            '* set ncarbonate = @ncarbonate ! Number of carbonate ions\n'
        ],
        [  # set OUTPUT (id set)
            'p1',
            '* set negid = CLA',
            '* set carbid = CO3\n'
        ],
        [  # set OUTPUT (id set)
            'p1',
            '* read sequence CLA',
            '* read sequence CO3 @ncarbonate !Generate CO3\n* generate CO3 first none last none setup\n'
        ]
    ]

    # write file "step2.2_ions.inp"
    patch_file(INP, list_replacements)

    # ===================================================== #
    # step2_solvator.inp

    INP = "step2_solvator.inp"

    list_replacements = [
        [  # define ncarbonate
            'a1',
            'stream step2.2_ions.prm\n',
            """
if ncarbonate .gt. 0 then
   read sequence @carbid @ncarbonate
   generate CO3 first none last none setup warn
endif
"""
        ],
        [  # calc nion
            's1',
            'calc nion = @npos + @nneg\n',
            'calc nion = @npos + @nneg + @ncarbonate\n'
        ],
        [  # set OUTPUT (number)
            'p1',
            '* set nneg = @nneg',
            '* set ncarbonate = @ncarbonate ! Number of carbonate ions\n'
        ],
        [  # set OUTPUT (id set)
            'p1',
            '* set negid = @negid',
            '* set carbid = @carbid\n'
        ],
        [  # set OUTPUT (define)
            'p1',
            '* if npos .gt. 0 then\n',
            """\
* if ncarbonate .gt. 0 then
*  read sequence @carbid @ncarbonate
*  generate CO3 first none last none setup warn
* endif
"""
        ]
    ]
    # write file
    patch_file(INP, list_replacements)

    # =====================================================#
    # step3_pbcsetup.inp
    INP = "step3_pbcsetup.inp"

    list_replacements = [
        [  # segid carbid
            's',
            'segid @posid .or. segid @negid',
            'segid @posid .or. segid @negid .or. segid @carbid'
        ],
        [  # set OUTPUT (set id)
            'p1',
            '* set negid = @negid',
            '* set carbid = @carbid\n'
        ]
    ]

    # write file
    patch_file(INP, list_replacements)
