#! usr/bin/env bash

# Run an interactive cloud workstation
dx run cloud_workstation \
    --allow-ssh \
    --brief \
    -y \
    -imax_session_length=4h \
    --priority high \
    --instance-type mem3_ssd1_v2_x8 \
    --name cloud_workstation