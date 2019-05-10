
    ima=double(imread('lena.png'));

    result=zeros(5,100);
  
    for le=1:5 %% 5 different noise levels
         level=le*10;
         seed   =  0;
         randn( 'state', seed );
        for i=1:100 %% 100 simulations for each noise level
            
            rima=ima+randn(size(ima))*level; %% add noise

            tic
            est=noiselevel(rima); %% noise level estimation
            res_t=toc;
            
            result(le,i)=est;
            R_t(le,i)=res_t;
 
        end
        
	disp(['Noise level : ' num2str(level) ', mean absolute error : ' num2str(mean(abs(result(le,:)-le*10))) ', maximum absolute error : ' num2str(max(abs(result(le,:)-le*10))) ]);

    end


