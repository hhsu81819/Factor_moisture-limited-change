clc
clear all
close all
addpath /homes/hhsu/Matlab_tool
addpath /homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis
load GSHHS_COAST220HL1
        GSHHS_COAST220HL2=GSHHS_COAST220HL1;
        GSHHS_COAST220HL1(:,1)=GSHHS_COAST220HL1(:,1)-180;
        GSHHS_COAST220HL2(:,1)=GSHHS_COAST220HL2(:,1)+180;
        MAP=cat(1,GSHHS_COAST220HL1,GSHHS_COAST220HL2);
        a=find(MAP(:,1)<-180);
        MAP(a,:)=[];
        a=find(MAP(:,1)>360);
        MAP(a,:)=[];
load /homes/hhsu/Colormap/weak_Blue_White_Red.mat
LAT=ncread('/project/cmip5/hhsu/AWI-ESM-1-1-LR/mrsos_AWI-ESM-1-1-LR_piControl_r1i1p1f1_regrided2x2_ng.nc','lat');
LON=ncread('/project/cmip5/hhsu/AWI-ESM-1-1-LR/mrsos_AWI-ESM-1-1-LR_piControl_r1i1p1f1_regrided2x2_ng.nc','lon');

alpha=0.05;

FileName={'AWI-ESM-1-1-LR','CanESM5','CMCC-ESM2','CNRM-CM6-1','CNRM-CM6-1-HR','INM-CM4-8','MIROC-ES2L','MRI-ESM2-0','CMCC-CM2-SR5','MPI-ESM-1-2-HAM','ICON-ESM-LR','IPSL-CM6A-LR','GFDL-CM4','NorESM2-MM','MIROC6'}


Var={'csm','wp','wepr','trpr','drpr','dpfc','dpfw','dcfp','dwfp','dcdp','dwdp'};
mpv1=nan(15,180,90);
mpv2=nan(15,180,90);
vote=nan(15,180,90);
        load 2d_rgb_colormap.mat

ccc=[142,1,82
197,27,125
222,119,174
241,182,218
253,224,239
230,245,208
184,225,134
127,188,65
77,146,33
39,100,25]/255;

ccc=flipud(ccc);

AA=permute(repmat(LAT,[1 180]),[2 1]);


for Model=1:15
        MODELNAME=char(FileName(Model))
	
	if exist(['/project/land/hhsu/04.CMIP6_trend/' MODELNAME '_Ptran.nc'],'file')
	SM=ncread(['/project/land/hhsu/04.CMIP6_trend/' MODELNAME '_SMcli_Pval.nc'],'SM_dif');
        SM(find(SM(:)>999))=nan;

	p1csm1=ncread(['/project/land/hhsu/04.CMIP6_trend/' MODELNAME '_Ptran.nc'],'p1csm1');
	p9csm1=ncread(['/project/land/hhsu/04.CMIP6_trend/' MODELNAME '_Ptran.nc'],'p9csm1');
	p1csm9=ncread(['/project/land/hhsu/04.CMIP6_trend/' MODELNAME '_Ptran.nc'],'p1csm9');
	p9csm9=ncread(['/project/land/hhsu/04.CMIP6_trend/' MODELNAME '_Ptran.nc'],'p9csm9');
	
	p9csm1(find(p9csm1(:)>999))=nan;
        p1csm1(find(p1csm1(:)>999))=nan;
        p1csm9(find(p1csm9(:)>999))=nan;
        p9csm9(find(p9csm9(:)>999))=nan;
	

	mpv1(Model,:,:)=squeeze(p9csm1-p1csm1)*100;
	mpv2(Model,:,:)=squeeze(p9csm9-p9csm1)*100; %
%
	for x=1:180
		for y=1:90
			if mpv1(Model,x,y)>-9999
				if y>33 & y<57
				DD=20*365;
				else
				DD=20*150;
				end
				kk=p9csm9(x,y);
				qq=p1csm1(x,y);
				[p, Q]= chi2test([kk*DD,DD-kk*DD;qq*DD,DD-qq*DD]);
				if p >0.05
				mpv1(Model,x,y)=nan;
				mpv2(Model,x,y)=nan;
				end
			end
		end
	end
%
	end
end
	
	for Model=1:15
	MODELNAME=char(FileName(Model))
	
	subplot(3,5,Model)
	pv1=squeeze(mpv1(Model,:,:));
	pv2=squeeze(mpv2(Model,:,:));

	sz=(pv1+pv2);
	SM=ncread(['/project/land/hhsu/04.CMIP6_trend/' MODELNAME '_SMcli_Pval.nc'],'SM_dif');
        SM(find(SM(:)>999))=nan;
	SM(find(SM(:)==0))=nan;
	SMo=SM;
	SMo(find(isnan(sz(:))))=nan;

	SM1=SMo;
        SM1(find(sz(:)>=0))=nan;
        SM1(find(SM1(:)>-9999))=1;
        SM1=SM1.*(abs(cos(AA*pi/2/90)));
        SM1=nansum(SM1(:));

        SM2=SMo;
        SM2(find(sz(:)<=0))=nan;
        SM2(find(SM2(:)>-9999))=1;
        SM2=SM2.*(abs(cos(AA*pi/2/90)));
        SM2=nansum(SM2(:));

        SM=SM;
        SM(find(SM(:)>-9999))=1;
        SM=SM.*(abs(cos(AA*pi/2/90)));
        ALL=nansum(SM(:));

        SM1=SM1/ALL*100;
        SM2=SM2/ALL*100;
        SM3=100-SM1-SM2;


	
	a1=sz(1:90,:);
        a2=sz(91:180,:);
        sz=cat(1,a2,a1);
	ssz=sz.^2;
	sz(find(ssz(:)<0.0000001))=nan;
	sz(90,:)=(sz(89,:)+sz(91,:))/2;

	sz(find(sz(:)==0))=nan;
        hold on
        pcolor(LON,LAT,sz'); shading flat
        caxis([-40 40])


   %     contourf(LON+1,LAT+1,sz',[-40:8:40],'LineColor','none')
%         pcolor(LON+1,LAT+1,sz'); shading flat
	hold on
%{
	for x=1:180
        for y=1:74

        if sz(x,y)< -40
        scatter(x*2,y*2-90,1,'Marker','d','MarkerEdgeColor', ccc(1,:),'MarkerFaceColor', ccc(1,:))
	end
        
	for a=1:8
        if sz(x,y)> -40+(a-1)*10 &  sz(x,y) < -40+a*10
        scatter(x*2,y*2-90,1,'Marker','d','MarkerEdgeColor',ccc(a+1,:),'MarkerFaceColor', ccc(a+1,:))
        end
        end

        if sz(x,y)> 40
        scatter(x*2,y*2-90,1,'Marker','d','MarkerEdgeColor',ccc(10,:),'MarkerFaceColor', ccc(10,:))
        end



        end
        end
%}
	


%        caxis([-50 50])
        colormap(ccc)
	hold on
	title({MODELNAME}) %,['\color[rgb]{0.1529    0.3922    0.0980}' num2str(round(SM1)) '%,\color[rgb]{0.5569    0.0039    0.3216}' num2str(round(SM2)) '%']})

        plot(MAP(:,1),MAP(:,2),'color',[0.7 0.7 0.7],'LineWidth',1.3);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',12)
        axis([50 340 -55 60])
end		
	ha=get(gcf,'children')
	 caxis([-40 40])
	colorbar('southoutside')

        ha=get(gcf,'children')
	ha
%{
        set(ha(1),'position',[.75 .03 .22 .02])
        set(ha(2),'position',[.01 .01 .24 .16])
        set(ha(3),'position',[.01 .22 .24 .16])
        set(ha(4),'position',[.01 .43 .24 .16])
        set(ha(5),'position',[.01 .64 .24 .16])
        set(ha(9),'position',[.25 .01 .24 .16])
        set(ha(7),'position',[.25 .22 .24 .16])
        set(ha(16),'position',[.25 .43 .24 .16])
        set(ha(10),'position',[.25 .64 .24 .16])
        set(ha(6),'position',[.49 .01 .24 .16])
        set(ha(11),'position',[.49 .22 .24 .16])
        set(ha(12),'position',[.49 .43 .24 .16])
        set(ha(13),'position',[.49 .64 .24 .16])
        set(ha(14),'position',[.73 .64 .24 .16])
        set(ha(15),'position',[.73  .22 .24 .16])
        set(ha(8),'position',[.73  .43 .24 .16])
%}
	set(ha(1),'position',[.3 .025 .4 .05])
        set(ha(9),'position',[.62 .82 .26 .14])
        set(ha(3),'position',[.36 .82 .26 .14])
        set(ha(4),'position',[.1 .82 .26 .14])
        set(ha(5),'position',[.1 .64 .26 .14])
        set(ha(6),'position',[.36 .64 .26 .14])
        set(ha(7),'position',[.62 .64 .26 .14])
        set(ha(15),'position',[.1 .46 .26 .14])
        set(ha(2),'position',[.36 .46 .26 .14])
        set(ha(10),'position',[.62 .46 .26 .14])
        set(ha(11),'position',[.1 .28 .26 .14])
        set(ha(12),'position',[.36 .28 .26 .14])
        set(ha(13),'position',[.62 .28 .26 .14])
        set(ha(16),'position',[.1  .1 .26 .14])
        set(ha(8),'position',[.36  .1 .26 .14])
        set(ha(14),'position',[.62  .1 .26 .14])


	set(gcf,'Units','centimeters','position',[1 1 40 30]);
        saveas(gcf,['./dPml_3x5_pcolor.png'])
        close all

%	set(gcf,'Units','centimeters','position',[1 1 45 20]);
  %      saveas(gcf,['./Composite_median_2dcolor_CSM_mos_eng_new.png'])
 %       close all




	
