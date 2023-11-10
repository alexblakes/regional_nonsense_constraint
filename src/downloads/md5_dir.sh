#!/usr/bin/env bash

# md5sum for all files in a directory

md5sum "${1}/*" > "${1}/${2}.md5"