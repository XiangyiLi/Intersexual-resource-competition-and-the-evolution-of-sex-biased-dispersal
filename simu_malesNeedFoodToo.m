% Introduction: 
% Note that some of the variables in the code have different names than in the main text, we clarify the meaning of each variable here:
% tmax: the maximum number of generation that the simulation runs
% nofsites: the number of sites
% phi: number of individuals in each site
% bm and bf: parameters that adjust the gradient between income/capital breeding
% mutrate: mutation rate
% resolution: time resolution of the dispersal period, is set to 20 in all simulations shown in the paper
% food: a parameter determines whether food competition is present between males and females in a patch
% malegreed: the relative food consumption of males compared to females, corresponding to "beta" in the main text of the paper


function []=simu_malesNeedFoodToo(tmax, nofsites, phi, bm, bf, mutrate, resolution, food, malegreed)
	
	function ind=ddists(p,n)
		ind=ones([1 n]);
		pp=cumsum(p/sum(p));
		rnd=rand([1 n]); 
		for i=1:length(pp)
			ind=ind+(rnd>pp(i));
        end
	end

    
	% N contains
	fd=[1 2]; % female dispersal alleles
	ft=[3 4]; % female timing alleles
	md=[5 6]; % male dispersal alleles
	mt=[7 8]; % male timing alleles
	sex=9;
	loc=10;
	condition=11;
	
	if length(food)>1
		DifferentPots=1;
	else
		DifferentPots=0;
	end


	DispCostf=(1-bf)*log(10); % dispersal cost parameter, right now the dependency on bm is there as this leads to 0.1 dispersal mortality for good condition individuals who've always been eating 1 unit of food per time unit
	DispCostm=(1-bm)*log(10);
	
	initpopsize = phi * nofsites;
	
	% initialize population
	N(1:initpopsize,fd)=0.5;
	N(1:initpopsize,ft)=resolution/2;
	N(1:initpopsize,md)=0.5;
	N(1:initpopsize,mt)=resolution/2;
	N(:,sex)=sign(randn(initpopsize,1)); % negative = males, positive = females
	N(:,loc)=unidrnd(nofsites,initpopsize,1);
	N(:,condition)=ones(initpopsize,1); % 1st generation - we assume everyone has same condition
	
	% initialize data collection
	allmeans=NaN(tmax,10);

	tmp=ones(phi,1)*(1:nofsites); tmp=tmp(:); Nnew=nan(size(tmp,1),condition); Nnew(:,loc)=tmp; Nnew(:,condition)=0; % this is the template for each new generation

	for t=1:tmax
		t
    
	    % collect genetic data
	    f=find(and(N(:,loc)>0, ~isnan(N(:,1))));
		for i=1:8
			allmeans(t,i)=mean(N(f,i));
		end
    
	    % new generation begins with breeding; who wins the site?
	    mother=nan(nofsites,1); father=nan(nofsites,1);
	    offspring=Nnew; % initialize new offspring for later

	    for k=1:nofsites
	        females=find(and(N(:,loc)==k, N(:,sex)>0)); % who is here: females
	        males=find(and(N(:,loc)==k, N(:,sex)<0)); % who is here: males
			
	        if and(length(females)>0, length(males)>0) % a breeding pair can be formed
	            mother=females(ddists(N(females,condition),1));
	            father=males(ddists(N(males,condition),1));
	            f=find(offspring(:,loc)==k); % whoever needs to be produced at this site has these individuals as mom and dad
	            offspring(f,fd(1))=N(mother,fd(unidrnd(2,phi,1)));
	            offspring(f,fd(2))=N(father,fd(unidrnd(2,phi,1)));
	            offspring(f,ft(1))=N(mother,ft(unidrnd(2,phi,1)));
	            offspring(f,ft(2))=N(father,ft(unidrnd(2,phi,1)));
	            offspring(f,md(1))=N(mother,md(unidrnd(2,phi,1)));
	            offspring(f,md(2))=N(father,md(unidrnd(2,phi,1)));
	            offspring(f,mt(1))=N(mother,mt(unidrnd(2,phi,1)));
	            offspring(f,mt(2))=N(father,mt(unidrnd(2,phi,1)));
	            offspring(f,sex)=sign(randn(phi,1));
	        end
	    end
		
	    % mutations
	    mut=rand(phi*nofsites,8)<mutrate;
	    offspring(:,fd)=(1-mut(:,fd)).*offspring(:,fd)+mut(:,fd).*rand(phi*nofsites,2);
	    offspring(:,md)=(1-mut(:,md)).*offspring(:,md)+mut(:,md).*rand(phi*nofsites,2);
		offspring(:,ft)=(1-mut(:,ft)).*offspring(:,ft)+mut(:,ft).*unidrnd(resolution,phi*nofsites,2);
	    offspring(:,mt)=(1-mut(:,mt)).*offspring(:,mt)+mut(:,mt).*unidrnd(resolution,phi*nofsites,2);
    
	    % discard the old population, replace with new
	    N=offspring;

	    % now check who will disperse
	    dispprop=mean(N(:,fd)')'.*(N(:,sex)>0)+mean(N(:,md)')'.*(N(:,sex)<0);
	    willdisp=rand(phi*nofsites,1)<dispprop;
	    % condition starts accumulating, and if someone wants to disperse it will do so
	    dispdata=zeros(resolution,4);
	    for j=1:resolution
	        % dispersal: any volunteers? Done separately for the 2 sexes simply to collect data on their condition differences at the time of dispersal
	        % female dispersers
	        ff=find(and(willdisp, N(:,sex)>0, round(mean(N(:,ft)')')==j)); dispdata(j,1)=length(ff); if dispdata(j,1)>0 dispdata(j,2)=mean(N(ff,condition)); end;
	        % male dispersers
	        fm=find(and(willdisp, N(:,sex)<0, round(mean(N(:,mt)')')==j)); dispdata(j,3)=length(fm); if dispdata(j,3)>0 dispdata(j,4)=mean(N(fm,condition)); end;
	        % all dispersers
	        if length(ff)>0
	            % if condition too poor to disperse, location becomes 0 (once pop reproduces, the program simply ignores individuals who don't have a positive integer location value)
	            doomedF=rand(length(ff),1)<exp(-DispCostf*N(ff,condition));
	            N(ff,loc)=unidrnd(nofsites,length(ff),1).*(1-doomedF); 
	        end
	        if length(fm)>0
	            % if condition too poor to disperse, location becomes 0 (once pop reproduces, the program simply ignores individuals who don't have a positive integer location value)
				doomedM=rand(length(fm),1)<exp(-DispCostm*N(fm,condition));
				N(fm,loc)=unidrnd(nofsites,length(fm),1).*(1-doomedM);% strictly speaking dispersal 'back' to own natal site isn't excluded, but this happens only with probability 1/nofsites
	        end
	        % then the condition accumulation, this'll be done for each site separately
	        for k=1:nofsites
                males = find(and(N(:,loc)==k, N(:,sex)<0));
                females = find(and(N(:,loc)==k, N(:,sex)>0));
				nfemales=length(females);
				nmales=length(males);
	            if DifferentPots==0 % males and females eat from the same pot
					if nmales>0 
						N(males,condition)=bm*N(males,condition)+food/(nfemales+malegreed*nmales); 
					end
					if nfemales>0 
						N(females,condition)=bf*N(females,condition)+food/(nfemales+malegreed*nmales); 
					end
	            else % males and females eat from different pots
					if nfemales>0 
						N(females,condition)=bf*N(females,condition)+food(1)/nfemales; 
					end
					if nmales>0 
						N(males,condition)=bm*N(males,condition)+food(2)/(nmales*malegreed); 
					end
	            end
	        end
	    end
	    % collect the data on mean condition at dispersal for females and males
	    allmeans(t,mt(end)+1)=sum(dispdata(:,1).*dispdata(:,2))/sum(dispdata(:,1)); % mean female condition at dispersal
	    allmeans(t,mt(end)+2)=sum(dispdata(:,3).*dispdata(:,4))/sum(dispdata(:,3)); % mean male condition at dispersal
	end
	csvwrite([strcat("MNFT-bf-", num2str(bf), "-bm-", num2str(bm),"-DifferentPots-", num2str(DifferentPots), "-malegreed-", num2str(malegreed), "-Output.csv")], allmeans);

end




















