#!/bin/bash

gcloud compute instances create st-alignment --project=velina-208320 --zone=us-central1-a --machine-type=e2-custom-2-8192 --network-interface=network-tier=PREMIUM,stack-type=IPV4_ONLY,subnet=managed --maintenance-policy=MIGRATE --provisioning-model=STANDARD --service-account=237653208032-compute@developer.gserviceaccount.com --scopes=https://www.googleapis.com/auth/devstorage.read_only,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/trace.append --enable-display-device --create-disk=auto-delete=yes,boot=yes,device-name=instance-4,image=projects/ubuntu-os-cloud/global/images/ubuntu-2004-focal-v20230302,mode=rw,size=40,type=projects/velina-208320/zones/us-central1-a/diskTypes/pd-balanced --no-shielded-secure-boot --shielded-vtpm --shielded-integrity-monitoring --labels=ec-src=vm_add-gcloud --reservation-affinity=any

echo "done"

