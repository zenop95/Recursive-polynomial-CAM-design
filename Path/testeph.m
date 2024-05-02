function testeph(efile,tfile)

% TESTEPH Test JPL ephemeris file reading and interpolation
%    TESTEPH(EFILE,TFILE) compares the results generated by reading
%    and interpolating the data in the binary ephemeris file EFILE
%    to the results from the ASCII test file TFILE.  The success
%    criterion is a difference of less than 1e-13.  A plot of the
%    difference is generated along with a message regarding test
%    success.
%
%    NOTE:  Some large libration values may be stored incorrectly
%    due to the number of significant digits involved.  Check any
%    errors in libration (NTARG=15) by comparing the 'Calculated
%    Result' column of the output against the reference value in
%    the ORIGINAL test file TFILE.

% Read the test data and initialize the ephemeris file

eph_global
testdat=read_testpo(tfile);
init_eph(efile)

ncases=size(testdat,1);

% Report how many cases were loaded and how many will be tested

disp(sprintf('\n%i cases loaded with Julian dates from %f to %f ...\n', ...
    ncases,min(testdat(:,1)),max(testdat(:,1))));

keep=find(testdat(:,1)>=SST(1) & testdat(:,1)<=SST(2));
testdat=testdat(keep,:);
ncases=size(testdat,1);

disp(sprintf('%i cases will be compared ...\n',size(keep,1))); 

% Pre-allocate variables for error and results

err=zeros(size(testdat,1),1);
result=zeros(size(testdat,1),1);

% For each case, interpolate the ephemeris and compare results

for i=1:ncases
    [r,v]=pleph(testdat(i,1),testdat(i,2),testdat(i,3),0);
    s=[r(:);v(:)];
    err(i)=s(testdat(i,4))-testdat(i,5);
    result(i)=s(testdat(i,4));
    if ~mod(i,200)
        disp(sprintf('%5i cases completed ...',i)) % Report progress
    end
end

% Set the tolerance and report success or detail failures

tol=1e-13;
ierr=find(abs(err)>tol);
if isempty(ierr)
    disp(sprintf('\nAll %i cases passed.  All errors are <= %5.0e\n',ncases,tol));
else
    disp(sprintf('\nThe following cases had errors > %5.0e\n',tol));
    disp(sprintf(' Case  Julian Date   ntarg ncent  dim    Reference Result    Calculated Result       Error'));
    disp(sprintf('----------------------------------------------------------------------------------------------'));
    disp(sprintf('%5i %12.2f %5i %5i %5i %20.13f %20.13f %15.4e\n', ...
        [ierr testdat(ierr,:) result(ierr) err(ierr)]'))
end

% Plot the results

figure
plot(testdat(:,1),err)
hold on
plot([testdat(1,1) testdat(size(testdat,1),1)],[tol tol],'r--')
plot([testdat(1,1) testdat(size(testdat,1),1)],[-tol -tol],'r--')
xlabel('Julian Date (days)')
ylabel('Error in Position or Velocity (AU,AU/day)')
title(sprintf('Ephemeris File: %s - Test File: %s',efile,tfile),'Interpreter','none');