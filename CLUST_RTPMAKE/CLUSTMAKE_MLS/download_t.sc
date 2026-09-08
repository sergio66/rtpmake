#!/bin/bash

GREP_OPTIONS=''

cookiejar=$(mktemp cookies.XXXXXXXXXX)
netrc=$(mktemp netrc.XXXXXXXXXX)
chmod 0600 "$cookiejar" "$netrc"
function finish {
  rm -rf "$cookiejar" "$netrc"
}

trap finish EXIT
WGETRC="$wgetrc"

prompt_credentials() {
    echo "Enter your Earthdata Login or other provider supplied credentials"
    read -p "Username (sergio66): " username
    username=${username:-sergio66}
    read -s -p "Password: " password
    echo "machine urs.earthdata.nasa.gov login $username password $password" >> $netrc
    echo
}

exit_with_error() {
    echo
    echo "Unable to Retrieve Data"
    echo
    echo $1
    echo
    echo "https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2025/MLS-Aura_L3MB-Temperature_v05-02-c01_2025.nc"
    echo
    exit 1
}

prompt_credentials
  detect_app_approval() {
    approved=`curl -s -b "$cookiejar" -c "$cookiejar" -L --max-redirs 5 --netrc-file "$netrc" https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2025/MLS-Aura_L3MB-Temperature_v05-02-c01_2025.nc -w '\n%{http_code}' | tail  -1`
    if [ "$approved" -ne "200" ] && [ "$approved" -ne "301" ] && [ "$approved" -ne "302" ]; then
        # User didn't approve the app. Direct users to approve the app in URS
        exit_with_error "Please ensure that you have authorized the remote application by visiting the link below "
    fi
}

setup_auth_curl() {
    # Firstly, check if it require URS authentication
    status=$(curl -s -z "$(date)" -w '\n%{http_code}' https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2025/MLS-Aura_L3MB-Temperature_v05-02-c01_2025.nc | tail -1)
    if [[ "$status" -ne "200" && "$status" -ne "304" ]]; then
        # URS authentication is required. Now further check if the application/remote service is approved.
        detect_app_approval
    fi
}

setup_auth_wget() {
    # The safest way to auth via curl is netrc. Note: there's no checking or feedback
    # if login is unsuccessful
    touch ~/.netrc
    chmod 0600 ~/.netrc
    credentials=$(grep 'machine urs.earthdata.nasa.gov' ~/.netrc)
    if [ -z "$credentials" ]; then
        cat "$netrc" >> ~/.netrc
    fi
}

fetch_urls() {
  if command -v curl >/dev/null 2>&1; then
      setup_auth_curl
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        curl -f -b "$cookiejar" -c "$cookiejar" -L --netrc-file "$netrc" -g -o $stripped_query_params -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  elif command -v wget >/dev/null 2>&1; then
      # We can't use wget to poke provider server to get info whether or not URS was integrated without download at least one of the files.
      echo
      echo "WARNING: Can't find curl, use wget instead."
      echo "WARNING: Script may not correctly identify Earthdata Login integrations."
      echo
      setup_auth_wget
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        wget --load-cookies "$cookiejar" --save-cookies "$cookiejar" --output-document $stripped_query_params --keep-session-cookies -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  else
      exit_with_error "Error: Could not find a command-line downloader.  Please install curl or wget"
  fi
}

fetch_urls <<'EDSCEOF'
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2025/MLS-Aura_L3MB-Temperature_v05-02-c01_2025.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2024/MLS-Aura_L3MB-Temperature_v05-02-c01_2024.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2023/MLS-Aura_L3MB-Temperature_v05-02-c01_2023.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2022/MLS-Aura_L3MB-Temperature_v05-02-c03_2022.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2021/MLS-Aura_L3MB-Temperature_v05-02-c02_2021.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2020/MLS-Aura_L3MB-Temperature_v05-02-c01_2020.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2019/MLS-Aura_L3MB-Temperature_v05-02-c01_2019.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2018/MLS-Aura_L3MB-Temperature_v05-02-c01_2018.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2017/MLS-Aura_L3MB-Temperature_v05-01-c06_2017.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2016/MLS-Aura_L3MB-Temperature_v05-01-c06_2016.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2015/MLS-Aura_L3MB-Temperature_v05-01-c06_2015.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2014/MLS-Aura_L3MB-Temperature_v05-01-c06_2014.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2013/MLS-Aura_L3MB-Temperature_v05-02-c01_2013.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2012/MLS-Aura_L3MB-Temperature_v05-02-c01_2012.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2011/MLS-Aura_L3MB-Temperature_v05-01-c06_2011.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2010/MLS-Aura_L3MB-Temperature_v05-01-c06_2010.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2009/MLS-Aura_L3MB-Temperature_v05-02-c01_2009.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2008/MLS-Aura_L3MB-Temperature_v05-01-c06_2008.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2007/MLS-Aura_L3MB-Temperature_v05-01-c06_2007.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2006/MLS-Aura_L3MB-Temperature_v05-02-c01_2006.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2005/MLS-Aura_L3MB-Temperature_v05-02-c01_2005.nc
https://data.gesdisc.earthdata.nasa.gov/data/Aura_MLS_Level3/ML3MBT.005/2004/MLS-Aura_L3MB-Temperature_v05-01-c06_2004.nc
EDSCEOF
