git config --global user.name Alex
git config --global user.email alex.blakes@gmail.com

eval "$(ssh-agent -s)"
ssh-add /opt/notebooks/id_ed25519

git clone git@github.com:alexblakes/ukb_constraint.git

cd ukb_constraint
eval "$(ssh-agent -s)"
ssh-add /opt/notebooks/id_ed25519
