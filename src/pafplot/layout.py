"""
#---------------------------------------------------------------- ParseIDs ----#
sub ParseIDs ($$$)
{
  my $aref = shift;
  my $rref = shift;
  my $qref = shift;

  my $align;

  foreach $align (@$aref)
  {

    my ($sR, $eR, $sQ, $eQ, $sim, $lenR, $lenQ, $idR, $idQ) = @$align;

    if ( !exists $rref->{$idR} )
    {
      $rref->{$idR} = [ $lenR - 1, $lenR, 1 ];
    }

    if ( !exists $qref->{$idQ} )
    {
      $qref->{$idQ} = [ $lenQ - 1, $lenQ, 1 ];
    }
  }
}
"""

def parse_ids(aref):
    rref = {}
    qref = {}
    for align in aref:
        if align.ref_id not in rref:
            rref[align.ref_id] = {
                "offset": None,
                "len": align.ref_length,
                "strand": 1
            }

        if align.query_id not in qref:
            qref[align.query_id] = {
                "offset": None,
                "len": align.query_length,
                "strand": 1
            }

    return rref, qref


"""
#--------------------------------------------------------------- LayoutIDs ----#
# For each reference and query sequence, find the set of alignments that
# produce the heaviest (both in non-redundant coverage and percent
# identity) alignment subset of each sequence using a modified version
# of the longest increasing subset algorithm. Let R be the union of all
# reference LIS subsets, and Q be the union of all query LIS
# subsets. Let S be the intersection of R and Q. Using this LIS subset,
# recursively span reference and query sequences by their smaller
# counterparts until all spanning sequences have been placed. The goal
# is to cluster all the "major" alignment information along the main
# diagonal for easy viewing and interpretation.
sub LayoutIDs ($$)
{
  my $rref = shift;
  my $qref = shift;

  my %rc;          # chains of qry seqs needed to span each ref
  my %qc;          # chains of ref seqs needed to span each qry
  #  {idR} -> [ placed, len, {idQ} -> [ \slope, \loR, \hiR, \loQ, \hiQ ] ]
  #  {idQ} -> [ placed, len, {idR} -> [ \slope, \loQ, \hiQ, \loR, \hiR ] ]

  my @rl;          # oo of ref seqs
  my @ql;          # oo of qry seqs
  #  [ [idR, slope] ]
  #  [ [idQ, slope] ]

  #-- get the filtered alignments
  open (MFILE, "<$OPT_Mfile")
    or die "ERROR: Could not open $OPT_Mfile, $!\n";

  my ($sR, $eR, $sQ, $eQ, $lenR, $lenQ, $idR, $idQ);
  my ($loR, $hiR, $loQ, $hiQ);
  my ($dR, $dQ, $slope);

  while ( <MFILE> ) {
    #--  >= 10 column match (Mashmap)
    if ( /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(([0-9]*[.])?[0-9]+).*$/ ) {

      $sR   = $8;   $eR   = $9;
      $lenR = $7;   $lenQ = $2;
      $idR  = $6;   $idQ  = $1;

      if ( $5 eq "+") {
        $sQ   = $3;   $eQ   = $4;
      } else {
        $sQ   = $4;   $eQ   = $3;
      }

      #-- skip it if not on include list
      if ( !exists $rref->{$idR} || !exists $qref->{$idQ} ) { next; }

      #-- get orientation of both alignments and alignment slope
      $dR = $sR < $eR ? 1 : -1;
      $dQ = $sQ < $eQ ? 1 : -1;
      $slope = $dR == $dQ ? 1 : -1;

      #-- get lo's and hi's
      $loR = $dR == 1 ? $sR : $eR;
      $hiR = $dR == 1 ? $eR : $sR;

      $loQ = $dQ == 1 ? $sQ : $eQ;
      $hiQ = $dQ == 1 ? $eQ : $sQ;

      if ($OPT_ONLY_USE_FATTEST)
      {
        #-- Check to see if there is another better alignment
        if (exists $qc{$idQ})
        {
          my ($oldR) = keys %{$qc{$idQ}[2]};
          my $val = $qc{$idQ}[2]{$oldR};

          if (${$val->[4]} - ${$val->[3]} > $hiR - $loR)
          {
            #-- Old alignment is better, skip this one
            next;
          }
          else
          {
            #-- This alignment is better, prune old alignment
            delete $rc{$oldR}[2]{$idQ};
            delete $qc{$idQ};
          }
        }
      }

      #-- initialize
      if ( !exists $rc{$idR} ) { $rc{$idR} = [ 0, $lenR, { } ]; }
      if ( !exists $qc{$idQ} ) { $qc{$idQ} = [ 0, $lenQ, { } ]; }

      #-- if no alignments for these two exist OR
      #-- this alignment is bigger than the current
      if ( !exists $rc{$idR}[2]{$idQ} || !exists $qc{$idQ}[2]{$idR} ||
        $hiR - $loR >
        ${$rc{$idR}[2]{$idQ}[2]} - ${$rc{$idR}[2]{$idQ}[1]} ) {

        #-- rc and qc reference these anonymous values
        my $aref = [ $slope, $loR, $hiR, $loQ, $hiQ ];

        #-- rc is ordered [ slope, loR, hiR, loQ, hiQ ]
        #-- qc is ordered [ slope, loQ, hiQ, loR, hiR ]
        $rc{$idR}[2]{$idQ}[0] = $qc{$idQ}[2]{$idR}[0] = \$aref->[0];
        $rc{$idR}[2]{$idQ}[1] = $qc{$idQ}[2]{$idR}[3] = \$aref->[1];
        $rc{$idR}[2]{$idQ}[2] = $qc{$idQ}[2]{$idR}[4] = \$aref->[2];
        $rc{$idR}[2]{$idQ}[3] = $qc{$idQ}[2]{$idR}[1] = \$aref->[3];
        $rc{$idR}[2]{$idQ}[4] = $qc{$idQ}[2]{$idR}[2] = \$aref->[4];
      }

      next;
    }

    #-- default
    die "ERROR: Could not parse $OPT_Mfile\n$_";
  }

  close (MFILE)
    or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";

  #-- recursively span sequences to generate the layout
  foreach $idR ( sort { $rc{$b}[1] <=> $rc{$a}[1] } keys %rc ) {
    SpanXwY ($idR, \%rc, \@rl, \%qc, \@ql);
  }

  #-- undefine the current offsets
  foreach $idR ( keys %{$rref} ) { undef $rref->{$idR}[0]; }
  foreach $idQ ( keys %{$qref} ) { undef $qref->{$idQ}[0]; }

  #-- redefine the offsets according to the new layout
  my $roff = 0;
  foreach my $r ( @rl ) {
    $idR = $r->[0];
    $rref->{$idR}[0] = $roff;
    $rref->{$idR}[2] = $r->[1];
    $roff += $rref->{$idR}[1] - 1;
  }
  #-- append the guys left out of the layout
  foreach $idR ( keys %{$rref} ) {
    if ( !defined $rref->{$idR}[0] ) {
      $rref->{$idR}[0] = $roff;
      $roff += $rref->{$idR}[1] - 1;
    }
  }

  #-- redefine the offsets according to the new layout
  my $qoff = 0;
  foreach my $q ( @ql ) {
    $idQ = $q->[0];
    $qref->{$idQ}[0] = $qoff;
    $qref->{$idQ}[2] = $q->[1];
    $qoff += $qref->{$idQ}[1] - 1;
  }
  #-- append the guys left out of the layout
  foreach $idQ ( keys %{$qref} ) {
    if ( !defined $qref->{$idQ}[0] ) {
      $qref->{$idQ}[0] = $qoff;
      $qoff += $qref->{$idQ}[1] - 1;
    }
  }
}
"""


def layout_ids(alignments, use_only_fattest=False):
    rc = dict()
    qc = dict()

    rrefs, qrefs = parse_ids(alignments)

    for line in alignments:
        idR = line.ref_id
        sR = line.ref_start
        eR = line.ref_end
        lenR = line.ref_length

        idQ = line.query_id
        sQ = line.query_start
        eQ = line.query_end
        lenQ = line.query_length

        dR = 1 if line.ref_strand == "+" else -1
        dQ = 1 if line.query_strand == "+" else -1
        slope = 1 if dR == dQ else -1

        loR = sR if dR == 1 else eR
        hiR = eR if dR == 1 else sR

        loQ = sQ if dQ == 1 else eQ
        hiQ = eQ if dQ == 1 else sQ

        # ? Use only fattest?
        if use_only_fattest and (idQ in qc):
            oldR = list(qc[idQ]["set"].keys())[0]
            val = qc[idQ]["set"][idR]
            if (val[4] - val[3]) > (hiR - loR):
                continue
            else:
                del rc[oldR]["set"][idQ]
                del qc[idQ]

        if idR not in rc:
            rc[idR] = {"is_placed": False, "len": lenR, "set": {}}

        if idQ not in qc:
            qc[idQ] = {"is_placed": False, "len": lenQ, "set": {}}

        if (
            (idQ not in rc[idR]["set"]) or
            (idR not in qc[idQ]["set"]) or
            ((hiR - loR) > (rc[idR]["set"][idQ][2] - rc[idR]["set"][idQ][1]))
        ):
            rc[idR]["set"][idQ] = (slope, loR, hiR, loQ, hiQ)
            qc[idQ]["set"][idR] = (slope, loQ, hiQ, loR, hiR)

    rl = []
    ql = []
    for idR in sorted(rc.keys(), key=lambda x: rc[x]["len"]):
        span_xwy(idR, rc, rl, qc, ql)

    for idR in rrefs.keys():
        rrefs[idR]["offset"] = None
    for idQ in qrefs.keys():
        qrefs[idQ]["offset"] = None

    roff = 0
    for idR, strandR  in rl:
        rrefs[idR]["offset"] = roff
        rrefs[idR]["strand"] = strandR

        roff += rrefs[idR]["len"] - 1

    for idR in rrefs.keys():
        if rrefs[idR]["offset"] is None:
            rrefs[idR]["offset"] = roff
            roff += rrefs[idR]["len"] - 1

    qoff = 0
    for idQ, strandQ in ql:
        qrefs[idQ]["offset"] = qoff
        qrefs[idQ]["strand"] = strandQ
        qoff += qrefs[idQ]["len"] - 1

    for idQ in qrefs.keys():
        if qrefs[idQ]["offset"] is None:
            qrefs[idQ]["offset"] = qoff
            roff += qrefs[idQ]["len"] - 1
    return rrefs, qrefs

"""
sub SpanXwY ($$$$$) {
    my $x   = shift;   # idX
    my $xcr = shift;   # xc ref
    my $xlr = shift;   # xl ref
    my $ycr = shift;   # yc ref
    my $ylr = shift;   # yl ref

    my @post;
    foreach my $y ( sort { ${$xcr->{$x}[2]{$a}[1]} <=> ${$xcr->{$x}[2]{$b}[1]} }
                    keys %{$xcr->{$x}[2]} ) {

        #-- skip if already placed (RECURSION BASE)
        if ( $ycr->{$y}[0] ) { next; }
        else { $ycr->{$y}[0] = 1; }

        #-- get len and slope info for y
        my $len = $ycr->{$y}[1];
        my $slope = ${$xcr->{$x}[2]{$y}[0]};

        #-- if we need to flip, reverse complement all y records
        if ( $slope == -1 ) {
            foreach my $xx ( keys %{$ycr->{$y}[2]} ) {
                ${$ycr->{$y}[2]{$xx}[0]} *= -1;

                my $loy = ${$ycr->{$y}[2]{$xx}[1]};
                my $hiy = ${$ycr->{$y}[2]{$xx}[2]};
                ${$ycr->{$y}[2]{$xx}[1]} = $len - $hiy + 1;
                ${$ycr->{$y}[2]{$xx}[2]} = $len - $loy + 1;
            }
        }

        #-- place y
        push @{$ylr}, [ $y, $slope ];

        #-- RECURSE if y > x, else save for later
        if ( $len > $xcr->{$x}[1] ) { SpanXwY ($y, $ycr, $ylr, $xcr, $xlr); }
        else { push @post, $y; }
    }

    #-- RECURSE for all y < x
    foreach my $y ( @post ) { SpanXwY ($y, $ycr, $ylr, $xcr, $xlr); }
}
"""


def span_xwy(x, xcr, ycr, xlr = None, ylr = None):
    if xlr is None:
        xlr = []

    if ylr is None:
        ylr = []

    post = []
    for y in sorted(xcr[x]["set"].keys(), key=lambda i: xcr[x]["set"][i][1]):
        if ycr[y]["is_placed"]:
            continue
        else:
            ycr[y]["is_placed"] = True

        len_ = ycr[y]["len"]
        slope = xcr[x]["set"][y][0]

        if slope == -1:
            for xx in ycr[y]["set"].keys():
                (sl, loy, hiy, lox, hix) = ycr[y]["set"][xx]
                sl *= -1

                tmp = len_ - hiy + 1
                hiy = len_ - loy + 1
                loy = tmp

                del tmp

                ycr[y]["set"][xx] = (sl, loy, hiy, lox, hix)

        ylr.append((y, slope))
        if len_ > xcr[x]["len"]:
            ylr, xlr = span_xwy(y, ycr, xcr, xlr = ylr, ylr = xlr)
        else:
            post.append(y)

    for y in post:
        ylr, xlr = span_xwy(y, ycr, xcr, xlr = ylr, ylr = xlr)

    return xlr, ylr
