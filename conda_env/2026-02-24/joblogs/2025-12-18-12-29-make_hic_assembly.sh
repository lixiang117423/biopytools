matlock bam2 juicer ../hic.filtered.bam out.links.mnd
sort -k2,2 -k6,6 out.links.mnd > out.sorted.links.mnd
python3 ~/software/3d-dna/utils/agp2assembly.py ../04.build/scaffolds.agp scaffolds.assembly
bash ~/software/3d-dna/visualize/run-asm-visualizer.sh -p false scaffolds.assembly out.sorted.links.mnd
