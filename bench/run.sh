nsyn=1000
cpr=1000

printf "%12s%12s%12s%12s%12s%12s%12s%12s\n" ranks run thread_tot advance comms walk enqueue exchange
for nranks in 1 2 4 8 16 32 64 128 256 1024 2048 4096 8192 16384
do
    param_file=params_$nranks.json
    out_file=out_$nranks

    cp params_base.json $param_file
    sed -i "s|NRANK|$nranks|g" $param_file
    sed -i "s|NSYN|$nsyn|g" $param_file
    sed -i "s|CPR|$cpr|g" $param_file
    drybench $param_file > $out_file

    run_t=`awk '/model-run/     {print $2}' $out_file`
    tot_t=`awk '/_p_ root/      {print $4}' $out_file`
    adv_t=`awk '/_p_   advance/ {print $4}' $out_file`
    com_t=`awk '/_p_   communication/ {print $4}' $out_file`
    wlk_t=`awk '/_p_     walkspikes/  {print $4}' $out_file`
    enq_t=`awk '/_p_     enqueue/     {print $4}'    $out_file`
    exg_t=`awk '/_p_     exchange/    {print $4}'   $out_file`

    printf "%12d%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n" $nranks $run_t $tot_t $adv_t $com_t $wlk_t $enq_t $exg_t
done


