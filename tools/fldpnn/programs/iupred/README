                       IUPred RELEASE NOTES
                       ====================

IUPred Version 1.0

by Zsuzsanna Dosztanyi


IUPred is supplied in source code form along with the required data files. The
program is written in ANSI C and the code should compile on any ANSI C compiler
e.g. the GNU C compiler. To be able to run outside the source directory, the
IUPred_PATH  environment variable has to be set to the location of the source
files. 

TO COMPILE: 

cc iupred.c -o iupred 

TO RUN IUPred:

iupred seqfile type
  
  where seqfile is the name of the sequence file
  
  type is any of the option of 

  	long
	short 
	glob
	
  for prediction of long disorder, short disorder ( e.g. missing residues in
  X-ray structures) or predicting globular domains. 


INPUT FILE: sequence_file in fasta format. One sequence per file.

EXAMPLE RUN: 

iupred P53_HUMAN.seq long


INTERPRETATION OF THE OUTPUT:

In the case of long and short types of disorder the output  gives the
likelihood of disorder for each residue, i.e. it is a value between 0 and 1,
and higher values indicate higher probability of disorder. Residues with values
above 0.5 can be regarded as disordered, and at this cutoff 5% of globular
proteins is expected to be predicted to disordered (false positives).
 
For the prediction type of globular domains it gives the number of globular
domains and list their start and end position in the sequence. This is followed
by the submitted sequence with residues of globular domains indicated by
uppercase letters. 


Please see the LICENSE file for the license terms for the software. It is
basically free for academic users, but a license fee applies to commercial
users. 

THE PUBLICATION OF RESEARCH USING IUPred MUST INCLUDE AN APPROPRIATE
CITATION TO THE METHOD:

The Pairwise Energy Content Estimated from Amino Acid Composition Discriminates 
between Folded and Intrinsically Unstructured Proteins
Zsuzsanna Dosztányi, Veronika Csizmók, Péter Tompa and István Simon
J. Mol. Biol. (2005) 347, 827-839.


SHORT SUMMARY OF THE METHOD

Intrinsically unstructured/disordered proteins have no single well-defined
tertiary structure in their native, functional state. Our server recognizes
such regions from the amino acid sequence based on the estimated pairwise
energy content. The underlying assumption is that globular proteins make a
large number of interresidue interactions, providing the stabilizing energy to
overcome the entropy loss during folding. In contrast, IUPs have special
sequences that do not have the capacity to form sufficient interresidue
interactions. Taking a set of globular proteins with known structure, we have
developed a simple formalism that allows the estimation of the pairwise
interaction energies of these proteins. It uses a quadratic expression in the
amino acid composition, which takes into account that the contribution of an
amino acid to order/disorder depends not only its own chemical type, but also
on its sequential environment, including its potential interaction partners.
Applying this calculation for IUP sequences, their estimated energies are
clearly shifted towards less favorable energies compared to globular proteins,
enabling the predicion of protein disorder on this ground. 
