# Track 2 methods — mechanism-constrained drug screen

PROBAND01. Not medical care, a diagnosis, or advice. Nothing here is a
recommendation to start, stop, or change a medicine.

**No market-approved medicine restores spindle-assembly-checkpoint
function.** The useful result is the filter that shows why the obvious
oncology hits are the wrong direction.

I am treating a biallelic *BUB1B* hypothesis as the starting point: a
stop-gain plus a last-exon missense (rank 1 in my Track 1 write-up).
The missense is still unproven and unphased. The drug argument stays
conditional on both copies being damaged.

Current care for mosaic variegated aneuploidy (MVA) is symptom
management and cancer surveillance. That does not change.


## Methods description

The official Track 2 methods form is voluntary except for the
generative-AI line. I filled the rest so the panel can see the method
without hunting through the narrative.

**Team.** bartholomewtj. Individual entry.

**Generative AI (required).** Used to develop genetic analysis tool, genome not run through any LLM model. Used for writing/structure of report only. (Claude Code via API)

**From mechanism to a candidate.** I wrote down what a valid drug
would have to do, then asked public drug-target databases what is
approved or late-phase on *BUB1B* and neighbouring checkpoint genes.
Hits were classed restore / compensate / worsen / trap by hand against
that test. Literature that sat outside the APIs (NMN / SIRT2, rapalogs,
HSCT, senolytics) was checked the same way. Nothing passed as an
approved restore.

**Automated or manual?** The gene-symbol queries are automated
(`scripts/screen.py`). Classing, the paediatric judgement, and the
literature rows are manual.

**Literature.** Full text for Suijkerbuijk 2010 (author manuscript),
North 2014 (PMC4194088), and Sieben 2020 (submitted version). Abstract
only for Malumbres 2024, Tang 2011, and Wang 2020. I did not invent
findings from papers I had not read in full.

**Public data.** Open Targets GraphQL, ChEMBL mechanisms, openFDA
labels, PubMed / PMC, Crossref. Snapshot tables in `scripts/out/`.

**Proprietary data.** None. The challenge genome was not used in this
screen.

**Mechanism.** Recessive loss of function at *BUB1B* (MVA syndrome 1).
The pathway is the spindle-assembly checkpoint. Downstream: chromosome
mis-segregation, mosaic aneuploidy, childhood cancer risk, and growth
restriction. Published biallelic genotypes are truncating plus missense
hypomorph. Putting BubR1 back restores the checkpoint in patient cells
(Suijkerbuijk 2010).

**Time.** Local Windows desktop. The API screen is minutes. No GPU.
No billed cloud compute. Initial imputation - ~8 hours via local PC.

**Abstract.** A harm-filtered gene-symbol screen for an approved
medicine that would restore spindle-assembly-checkpoint function in
biallelic *BUB1B*. There isn't one. The oncology hits on these genes
weaken the checkpoint or need a working checkpoint to act. That is the
finding. Strength: the valid-drug test is written down before the
databases are queried, so wrong-direction compounds cannot sneak in as
candidates. Limitation: no patient-cell assay in this hackathon; the
Track 1 missense is unproven and unphased; the argument is conditional
on biallelic *BUB1B*.


## What a valid drug would have to do

MVA syndrome 1 is recessive loss of function. Published biallelic
genotypes pair a truncating allele with a missense allele. Patient
cells have low BubR1 (the *BUB1B* protein), a weak mitotic checkpoint,
and chromosome alignment defects. Putting BubR1 back restores the
checkpoint (Suijkerbuijk et al., 2010,
doi:10.1158/0008-5472.CAN-09-4319, author manuscript via PMC).

So a candidate has to:

1. **Raise checkpoint strength** in patient-derived cells — not merely
   kill aneuploid cells.
2. **Not make mitotic fidelity worse.** Childhood cancer risk is part of
   MVA1. Checkpoint-killing oncology drugs push the wrong way.
3. **Not trade growth for a theoretical mitotic delay.** Growth
   restriction is already in the syndrome. That rules out mTOR
   inhibitors even though mouse papers mention mTOR.

Kinase-domain missenses in published MVA1 act as hypomorphs by lowering
BubR1 protein 2–6-fold; HSP90 inhibition drops them further
(geldanamycin, Suijkerbuijk Fig. 3C). The Track 1 last-exon missense
sits in that C-terminal kinase neighbourhood, ten residues from the
established hypomorph p.Leu1012Pro, and well away from the
kinase-domain WT-like allele p.Gln921His. Neighbourhood is not an
assay. Until abundance or replacement data exist, it remains a
plausible, unproven hypomorph. Mouse BubR1 L1002P models human
L1012P, not this allele.


## Screen

I queried public APIs using gene symbols only — no variant coordinates,
no phenotype text, no VCF.

Script: `scripts/screen.py`. Snapshot: `scripts/out/`.

- Open Targets GraphQL (`drugAndClinicalCandidates`, tractability,
  chemical probes)
- ChEMBL mechanisms
- openFDA labels for unique approved names (capped at 40)

Genes: *BUB1B* and SAC neighbours (*BUB1*, *BUB3*, *MAD2L1*, *MAD1L1*,
*CDC20*); other MVA genes (*CEP57*, *TRIP13*, *MAD2L1BP*); then the
nodes that look tempting and are usually the wrong direction
(*AURKB*, *PLK1*, *TTK*, *TUBB*, *HSP90AA1*, *MTOR*, *SIRT2*).

Each hit is classed restore / compensate / worsen / trap.


## Result

| Gene | Clinical candidates (Open Targets) | Approved on this field |
|---|---:|---|
| *BUB1B* | 0 | 0 |
| *BUB1*, *BUB3*, *MAD2L1*, *MAD1L1*, *CDC20* | 0 | 0 |
| *CEP57*, *TRIP13*, *MAD2L1BP* | 0 | 0 |
| *SIRT2* | 0 | 0 |
| *AURKB* | 18 | 0 |
| *PLK1* | 8 | 0 |
| *TTK* | 2 | 0 |
| *HSP90AA1* | 10 | 0 |
| *MTOR* | 26 | 2 (perhexiline, ridaforolimus) |
| *TUBB* | 88 | 25 microtubule poisons / ADCs |

Full rows: `scripts/out/genes.tsv` and `hits.tsv`.

The same empty pattern showed up on the other MVA genes. The screen
only needs gene symbols, so it can be rerun for another case without
touching a genome.


## Reject table (this is the finding)

| Class | Examples | Why not |
|---|---|---|
| Worsen — SAC / mitotic kinases | Aurora B (barasertib, tozasertib, ibcasertib), PLK1 (volasertib, onvansertib), MPS1/TTK (BAY-1217389). None approved. | They weaken the checkpoint. MVA1 already has a weak checkpoint. |
| Worsen — HSP90 | Tanespimycin, ganetespib, retaspimycin. None approved. | Suijkerbuijk: HSP90 inhibition further drops kinase-domain missense BubR1 (geldanamycin, Fig. 3C). The opposite of restore. |
| Compensate-toxic — microtubules | Paclitaxel, docetaxel, vincristine, vinblastine, vinorelbine, colchicine, eribulin, several vedotin ADCs. Approved. | These drugs need a working SAC to hold mitotic arrest. MVA patient cells fail that arrest on nocodazole / STLC. They are also childhood cytotoxics on a cancer-predisposition background. |
| Trap — mTOR | Rapalogs (sirolimus, everolimus; used in children with TSC). Open Targets' MTOR list on this pull did not return those two, but they are approved for other indications. Ridaforolimus was listed. | Sieben et al. 2019/2020 (doi:10.1172/JCI126863, submitted version) linked monoallelic BubR1 mouse sarcopenia to mTORC1 hyperactivity and cited aging-rapamycin papers. They did not treat those mice. Rapalogs impair growth. MVA already includes growth restriction. Not a SAC restore. |

I am not proposing any of these as leads.


## Candidate table

| Rank | Drug | Class | Confidence |
|---|---|---|---|
| — | *none* | restore, approved | — |

Track 2 asked for existing market-approved medicines. I do not have
one. What I looked at beyond the gene-symbol screen, with a DOI and a
why / why-not, is in
[bartholomewtj_track2_considered.md](bartholomewtj_track2_considered.md).
None of those rows is a candidate.

One literature near-miss, **not** a candidate: North et al., 2014
(doi:10.15252/embj.201386907, PMC4194088). Full text.

- Mechanism: CBP acetylates BubR1 at K668 and sends it to the
  proteasome; SIRT2 deacetylates K668 and stabilises the protein.
- NMN (500 mg/kg/day IP × 7, Sigma reagent) raised NAD+ and restored
  BubR1 protein in **aged wild-type** mouse testes, SIRT2-dependent in
  MEFs. That experiment is not BubR1-hypomorphic mice and not patient
  cells.
- SIRT2 overexpression on BubR1H/H mice raised median lifespan 58%
  (males 123%; females unchanged) and partly rescued heart size /
  J-point. Lordokyphosis, cataracts, and small size did not reverse.
  Genetic overexpression is not a drug.
- Nicotinamide and sirtinol **drop** BubR1. NAM is not interchangeable
  with NMN.

NMN is not a registered medicine. *SIRT2* had no approved activator on
this screen. I am not treating it as a Track 2 candidate.

I could not get the full text of Malumbres & Villarroya-Beltri 2024
(doi:10.1038/s41576-024-00762-6, *Nat Rev Genet*). The abstract
(PMID 39169218) names no drug. The reference list includes
Suijkerbuijk, the JCI mouse paper, an HSCT case report, and CIN /
chemo-resistance papers — not a restore-class approved medicine.


## Next experiment

I did not treat cells.

The experiment that would matter: patient-derived fibroblasts or
lymphoblasts, mitotic delay on nocodazole or STLC, and BubR1
abundance, with and without a future restore-class compound.
Suijkerbuijk already showed that putting BubR1 back works. A drug has
to copy that readout.

Cancer surveillance stays the clinical standard. This packet does not
change that.


## Scalability

The screen is gene symbols in, classified table out. No genome. The
same run is empty on *CEP57*, *TRIP13*, and *MAD2L1BP*. That is the
reusable piece for other N-of-1 SAC disorders: lock the mechanism,
dump the wrong-direction mitosis drugs, and only then look for an
approved restore.

An empty candidate table is still useful for this child. The oncology
list on these genes is inhibitors and microtubule poisons. Those would
make a weak checkpoint weaker, or they need a working checkpoint to
do their job. Writing that down is the thing that could stop a
repurposing attempt going the wrong way.


## Acknowledgement

This work was made possible through the Hackathon, organized by Sage
Bionetworks in partnership with the MVA Society, Hugging Face, and
BEACON (The Benchmarking, Evaluation, and Assessment Consortium for
Science), with prize sponsorship from AWS and Anthropic. We are deeply
grateful to the child and their family who generously contributed their
data and their story to advance research into this rare disease. We
acknowledge their trust in making this Hackathon possible.
