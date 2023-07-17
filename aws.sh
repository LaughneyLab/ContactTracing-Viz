#!/usr/bin/env sh

# Make sure you have the elastic beanstalk cli installed (pipx install awsebcli)
# Make sure ssh is set up with eb ssh --setup
# Also make sure you have the AWS CLI installed https://aws.amazon.com/cli/
# The application on EBS is configured to use the following server hardware:
# 1. t3a.large (2 vCPU, 8GB RAM)
# 2. t3.large (2 vCPU, 8GB RAM)

# When datasets are updated and figures are precompiled, you must update the Amazon EFS
# Use https://us-east-2.console.aws.amazon.com/transfer#/ to create and instance so that you can transfer the data to and from
# https://docs.aws.amazon.com/transfer/latest/userguide/transfer-file.html#openssh
# Enter this directory, run the following command:
# sftp -i s-xxxx yyyy@s-xxxxb.server.transfer.us-east-2.amazonaws.com
# Then run: put -r data

# Next, you must update the EFS
# The .ebextensions script should mount it to my_efs
# More info here: https://github.com/oreaba/scripts/tree/main/mount_efs_beanstalk
# Additionally there is a ebextensions script to allow for more time for the docker container to build

# Since Docker images require a lot of compute to build, we build them locally using the Dockerfile and
# Upload them to our public ECR (Elastic Container Registry)

# For first time setup, run:
aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws/d6q7u0s0

# Build locally first
docker build --platform linux/amd64 --compress --rm -t ct-site .

# To run locally: docker run -v "$(pwd)/data":/app/data -p 8000:8000 ct-site

# Tag the build
docker tag ct-site:latest public.ecr.aws/d6q7u0s0/ct-site:latest

# Push the build to ECR
docker push public.ecr.aws/d6q7u0s0/ct-site:latest

# Deploy the build to EBS
eb deploy --timeout 30

# Clean all local docker data
docker system prune -a

# We can follow the logs with:
# eb ssh -c "tail -f /var/log/eb-engine.log"
# Or after deployment attempt, download all log files and view them locally

eb open