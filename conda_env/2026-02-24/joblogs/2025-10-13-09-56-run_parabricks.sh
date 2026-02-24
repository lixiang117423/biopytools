bwa index  /share/org/YZWL/yzwl_lixg/tmp/parabricks/genome.fa


biopytools parabricks \
  -i /share/org/YZWL/yzwl_lixg/tmp/parabricks/data \
  -o /share/org/YZWL/yzwl_lixg/tmp/parabricks/mapping \
  -r /share/org/YZWL/yzwl_lixg/tmp/parabricks/genome.fa \
  -t 64
