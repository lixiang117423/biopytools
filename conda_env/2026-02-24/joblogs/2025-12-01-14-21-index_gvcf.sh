# 方案一：使用for循环（简单直接）
for file in /share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint/vcf/*.gz; do
    tabix -@ -p vcf "$file"
done
