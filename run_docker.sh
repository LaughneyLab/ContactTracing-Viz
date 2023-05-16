#!/usr/bin/env sh

# Add -d to daemonize
docker run --rm --name ct_viz -p 8000:8000 -v data:/app/data ct-site