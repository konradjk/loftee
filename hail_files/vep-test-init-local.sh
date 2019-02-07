#!/bin/bash

VEP_LOC=/Users/konradk/vep
PERL_LIB_DIR=/opt/local/lib/perl5/5.26/

# Copy VEP
mkdir -p $VEP_LOC/homo_sapiens 
gsutil -m cp -r gs://konradk/vep/loftee $VEP_LOC
gsutil -m cp -r gs://hail-common/vep/vep/ensembl-tools-release-85 $VEP_LOC
gsutil -m cp -r gs://hail-common/vep/vep/loftee_data $VEP_LOC
gsutil -m cp -r gs://hail-common/vep/vep/Plugins $VEP_LOC
gsutil -m cp -r gs://hail-common/vep/vep/homo_sapiens/85_GRCh37 $VEP_LOC/homo_sapiens/
gsutil cp gs://hail-common/vep/vep/vep85-gcloud.json $VEP_LOC/vep85-gcloud.json
gsutil cp gs://hail-common/vep/vep/vep85-gcloud.properties $VEP_LOC/vep-gcloud.properties
gsutil cp gs://konradk/vep/phylocsf_gerp.sql $VEP_LOC/loftee_data/phylocsf.sql

#Create symlink to vep
ln -s $VEP_LOC/ensembl-tools-release-85/scripts/variant_effect_predictor $VEP_LOC

#Give perms
chmod -R 777 $VEP_LOC

# Copy perl JSON module
gsutil -m cp -r gs://hail-common/vep/perl-JSON/* $PERL_LIB_DIR

#Copy perl DBD::SQLite module
gsutil -m cp -r gs://hail-common/vep/perl-SQLITE/* $PERL_LIB_DIR

gsutil -m cp -r gs://konradk/vep/perl-libs/* $PERL_LIB_DIR

#Run VEP on the 1-variant VCF to create fasta.index file -- caution do not make fasta.index file writeable afterwards!
gsutil cp gs://konradk/vep/1var.vcf $VEP_LOC
gsutil cp gs://hail-common/vep/vep/run_hail_vep85_vcf.sh $VEP_LOC
chmod a+rx $VEP_LOC/run_hail_vep85_vcf.sh

$VEP_LOC/run_hail_vep85_vcf.sh $VEP_LOC/1var.vcf
