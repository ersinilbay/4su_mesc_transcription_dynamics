"""
Curated gene lists used for figure overlays and interpretation.

Important:
These grouped genes are not arbitrary marker sets. They were organized post hoc
from downstream candidate-gene output, including NASC-seq2-derived results,
to make residual-Fano and related figures more interpretable.
"""

GENES_REP1 = [
    "2410141K09Rik","Abcf2","Api5","Aplp2","Arpc1a","B230219D22Rik","Ccdc59","Cdca8","Cdk9",
    "Cry1","Csnk1a1","Ddx47","Dnajc21","Eda2r","Eps15l1","Faf1","Fam168b","Farsa","Fnta","Gar1",
    "Gatad2a","Gid8","Gm42418","Gpi1","Hmmr","Hspd1","Isg20l2","Jade1","Kat6a","Keap1","Lpin1",
    "Lsm14a","Malat1","Mat2b","Mob4","mt-Nd1","mt-Nd5","Nanog","Nom1","Nus1","Paxip1","Pbk",
    "Ptbp1","Rad18","Ranbp10","Rbm25","Rfx7","Rpl3","Scpep1","Sin3b","Smn1","Spice1","Tdgf1",
    "Thumpd3","Top2a","Tra2a","Trim33","Uso1","Wbp4","Xiap","Xrcc5","Zfp106","Zfp42"
]

GENES_REP2 = [
    "Cct7","Cenpj","Copb1","Cyb5r3","Dag1","Ddx47","Dnajc5","Edrf1","Eif3a","Fam193a",
    "Gm12346","Gm26624","Gm42418","Golga7","H3f3b","Hnrnph1","Hspa14","Incenp","Jade1",
    "Map3k7","Map4","Mrps31","Mybl2","Ncor1","Nup205","Pnn","Prpf3","Prpf8","Rbbp6",
    "Rhebl1","Rpl23","Rtf1","Scpep1","Slc38a1","Smc4","Tnks","Tpt1","Txndc9","Uso1",
    "Wtap","Zfp654","Zic3"
]

SELECTED_LABELS = {
    "Pluripotent_rep1": ["Nanog","Zfp42","Tdgf1","Rpl3","Farsa","Gpi1","Cdk9"],
    "Pluripotent_rep2": ["Zic3","Mybl2","Ncor1","Eif3a","Rpl23","H3f3b"],
}

REP1_PLURI_GROUPS = {
    "Pluripotency / early-embryo signaling": [
        "Zfp42","Dppa5a","Lefty2","L1td1","Tdgf1"
    ],
    "Chromatin / genome regulation": [
        "Set","Top2a","Smc4","H3f3b","Smarca5","Hmgb2","Nasp","Nap1l1",
        "Mtf2","Pds5a"
    ],
    "RNA biology": [
        "Dqx1","Hnrnpa2b1","Hnrnpc","G3bp2","Eif2s2","Eif3a","Eif5a",
        "Eef1a1","Eef1b2","Eef2","Nop56","Nop58","Rbm25","Neat1",
        "Mrps31","Rpl23"
    ],
    "Other regulators": [
        "Ankfy1","Arl6ip1","Cbfb","Cep95","Ctsl","Dennd5b","Fkbp3","Hmmr",
        "Macrod2","Malat1","Mdm2","Nedd4","Nek1","Npl","Nucks1","Pdgfrl",
        "Peg10","Psma7","Ran","Samd4b","Sgsm3","mt-Nd5"
    ],
}

REP1_GROUP_COLORS = {
    "Pluripotency / early-embryo signaling": "#2ca02c",
    "Chromatin / genome regulation": "#ff7f0e",
    "RNA biology": "#9467bd",
    "Other regulators": "crimson",
}

REP1_GROUP_MARKERS = {
    "Pluripotency / early-embryo signaling": "o",
    "Chromatin / genome regulation": "s",
    "RNA biology": "D",
    "Other regulators": "o",
}

REP2_PLURI_GROUPS = {
    "Pluripotency / early-embryo signaling": [
        "Zfp42","Dppa5a","Lefty2","L1td1","Tdgf1"
    ],
    "Chromatin / genome regulation": [
        "Top2a","Smc4","H3f3b","Mtf2","Nap1l1","Hmgb2",
        "Cenpf","Pcna","Ube2c","Mdm2"
    ],
    "RNA biology": [
        "Dqx1","Hnrnpa2b1","Hnrnpc","Hnrnpm","Hnrnpab","G3bp2",
        "Eif3a","Eif2s2","Eif5a","Eef1a1","Eef1b2","Eef2",
        "Nop56","Nop58","Rbm25","Neat1","Rsl1d1","Pcbp2","Cct6a","Rpl23"
    ],
    "Other regulators": [
        "Fkbp3","Psma7","Hsp90b1","Hspd1","Nedd4","Nucks1","Npl","Ran",
        "Pdgfrl","Peg10","mt-Nd5",
        "Ldha","Pkm","Atp5b","Cox4i1","Slc25a5","Gpx1","Prdx1","Glrx2","Mt1","Mt2",
        "Vim","Tpm1","Tpm3","Tuba1b","Add1","Sparc","Fn1","Cald1",
        "Apoe","Fabp3","Nefl","Etfdh","Ddit4","Pmaip1","Calm1"
    ],
}

REP2_GROUP_COLORS = REP1_GROUP_COLORS
REP2_GROUP_MARKERS = REP1_GROUP_MARKERS

PLURI_REP1_BOTTOMRIGHT = [
    "2410141K09Rik","Abcf2","Api5","Aplp2","Arpc1a","B230219D22Rik","Ccdc59","Cdca8","Cdk9","Ddx47","Eps15l1","Faf1",
    "Fam168b","Farsa","Fnta","Gar1","Gatad2a","Gid8","Gpi1","Isg20l2","Keap1","Lpin1","Lsm14a","Mat2b",
    "Mob4","Nom1","Nus1","Paxip1","Ptbp1","Ranbp10","Rfx7","Rpl3","Scpep1","Sin3b","Smn1","Thumpd3",
    "Uso1","Wbp4","Xiap"
]
PLURI_REP1_TOPLEFT = [
    "Cry1","Csnk1a1","Dnajc21","Eda2r","Gm42418","Hmmr","Hspd1","Jade1","Kat6a","Malat1","mt-Nd1","mt-Nd5",
    "Nanog","Rad18","Rbm25","Spice1","Tdgf1","Top2a","Tra2a","Trim33","Xrcc5","Zfp106","Zfp42"
]
PLURI_REP2_BOTTOMRIGHT = [
    "Cct7","Cenpj","Copb1","Cyb5r3","Dag1","Ddx47","Dnajc5","Edrf1","Fam193a","Gm12346","Golga7","Hspa14",
    "Map3k7","Map4","Mybl2","Nup205","Prpf3","Prpf8","Rhebl1","Rpl23","Rtf1","Scpep1","Slc38a1","Tnks",
    "Tpt1","Txndc9","Uso1","Wtap","Zfp654","Zic3"
]
PLURI_REP2_TOPLEFT = [
    "Eif3a","Gm26624","Gm42418","H3f3b","Hnrnph1","Incenp","Jade1","Mrps31","Ncor1","Pnn","Rbbp6","Smc4"
]

GROUP_COLORS = {"highFreq_lowSize": "crimson", "lowFreq_highSize": "purple"}
GROUP_MARKERS = {"highFreq_lowSize": "o", "lowFreq_highSize": "D"}

DIST_GENES_PLURI = ["Tdgf1", "Farsa"]

PRIORITY_2C = ["Zscan4","Zscan4c","Zscan4d","Zscan4f","Zscan4-ps1","Zscan4-ps2","Zscan4-ps3"]
PRIORITY_LENIENT = dict(min_mean=1e-3, max_mean=5.0, min_fold=1.7)