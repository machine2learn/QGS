for dir in */ ; do
  cd "$dir"
  res=`python3 test.py 2>&1`
  retval=$?
  if [ "$retval" -eq "0" ]; then
    echo "Test $dir successful."
  else
    echo "Test $dir FAILED."
    echo "$res"
  fi
  cd ..
done
