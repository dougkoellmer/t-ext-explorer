cd ../

rm -rf ./front_end_build
mkdir ./front_end_build

cp index.html ./front_end_build

cp -r ./src ./front_end_build
cp -r ./lib ./front_end_build
cp -r ./img ./front_end_build

zip -r ./front_end_build.zip ./front_end_build