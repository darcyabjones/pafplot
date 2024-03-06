def layout_ids(rrefs, qrefs, alignments, use_only_fattest=False):
    rc = dict()
    qc = dict()

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
            oldR = list(qc[idQ][2].keys())[0]
            val = qc[idQ][2][idR]
            if (val[4] - val[3]) > (hiR - loR):
                continue
            else:
                del rc[oldR][2][idQ]
                del qc[idQ]

        if idR not in rc:
            rc[idR] = (0, lenR, {})

        if idQ not in qc:
            qc[idQ] = (0, lenQ, {})

        if (
            (idQ not in rc[idR][2]) or
            (idR not in qc[idQ][2]) or
            ((hiR - loR) > (rc[idR][2][idQ][2] - rc[idR][2][idQ][1]))
        ):
            rc[idR][2][idQ] = (slope, loR, hiR, loQ, hiQ)
            qc[idQ][2][idR] = (slope, loQ, hiQ, loR, hiR)


    # Not finished
    return


def span_xwy():
    return
