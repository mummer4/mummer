/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

while(True)
{
  if(stoptraversal != NULL && stoptraversal(stopinfo))
  {
    return 0;
  }
  if(currentnode.toleaf)
  {
    DEBUG1(3,"visit leaf %lu ",
              (Showuint) LEAFADDR2NUM(stree,currentnode.address));
    DEBUG1(3,"below %lu\n",(Showuint) BRADDR2NUM(stree,startnode->address));
    if(processleaf(LEAFADDR2NUM(stree,currentnode.address),lcpnode,info) != 0)
    {
      return -1;
    }
    brotherval = LEAFBROTHERVAL(*(currentnode.address));
    if(NILPTR(brotherval))
    {
      readyforpop = True;
      currentnode.toleaf = False;
    } else
    {
      SETCURRENT(brotherval);     // current comes from brother
      lcpnode = stack.spaceBref[stack.nextfreeBref-1];
    }
  } else
  {
    if(readyforpop)
    {
      if(stack.nextfreeBref == UintConst(1))
      {
        break;
      }
      (stack.nextfreeBref)--;
      DEBUG1(3,"#pop[%lu]=",(Showuint) stack.nextfreeBref);
      DEBUG1(3,"%lu\n",
             (Showuint) BRADDR2NUM(stree,stack.spaceBref[stack.nextfreeBref]));
      PROCESSBRANCH2(stack.spaceBref[stack.nextfreeBref],info);
      brotherval = GETBROTHER(stack.spaceBref[stack.nextfreeBref]);
      if(!NILPTR(brotherval))
      {
        SETCURRENT(brotherval);    // current comes from brother
        lcpnode = stack.spaceBref[stack.nextfreeBref-1];
        readyforpop = False;
      }
    } else
    {
      DEBUG1(3,"#process1 %lu\n",
               (Showuint) BRADDR2NUM(stree,currentnode.address));
      PROCESSBRANCH1(currentnode.address,info);
      if(godown)
      {
        STOREINARRAY(&stack,Bref,128,currentnode.address);
        DEBUG1(3,"#push[%lu]=",(Showuint) (stack.nextfreeBref-1));
        DEBUG1(3,"%lu\n",(Showuint) BRADDR2NUM(stree,currentnode.address));
        child = GETCHILD(currentnode.address);
        SETCURRENT(child);    // current comes from child
      } else
      {
        brotherval = GETBROTHER(currentnode.address);
        if(NILPTR(brotherval))
        {
          readyforpop = True;
        } else
        {
          SETCURRENT(brotherval);    // current comes brother
        }
      }
    }
  }
}
