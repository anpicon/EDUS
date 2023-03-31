fp_Loss.close(); fp_E.close();
if(iCurrent)
{
    fp_J1.close(); fp_J2.close(); fp_J3.close(); fp_Jintra.close();
}
if(iTAbs)
{
    fp_TAbs.close();
}

time_t time_ = time(0);
tm* now_ = localtime(&time_);

printf("Calculation completed day %4d/%02d/%02d at %02d.%02d.%02d \n", (now_->tm_year + 1900), (now_->tm_mon + 1), (now_->tm_mday), (now_->tm_hour), (now_->tm_min), (now_->tm_sec));

int deltad =  (now_->tm_mday) - iday;
int deltah =  (now_->tm_hour) - ihour;

int deltam =  (now_->tm_min)  - imin;
int deltas =  (now_->tm_sec)  - isec;


if (deltas < 0)
{
    deltas += 60;
    deltam -= 1;
}
if (deltam < 0)
{
    deltam += 60;
    deltah -= 1;
}
if (deltah < 0)
{
    deltah += 24;
    deltad -= 1;
}


printf("Elapsed time:  %02d days, %02d hours, %02d minutes and %02d seconds.\n", deltad, deltah, deltam, deltas);

clock2=clock();
printf("Time of computation: %8.4f s\n", (double) (clock2-clock1)/CLOCKS_PER_SEC);

return 0;
